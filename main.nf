// Needed to publish results (NF 24.10.x preview feature)
nextflow.preview.output = true

// import the subworkflows
include { LOCAL_GET_CONTIG_INFO           } from './modules/local/get_contig_info/main'
include { LOCAL_COMPILE_EIGENS            } from './modules/local/compile_eigens/main'
include { VALIDATE_REPLICATION_COMPLETE   } from './modules/local/validate_replication_complete/main'
include { DB_CLEAN_REPLAY_SLOTS           } from './modules/local/db_clean_replay_slots/main'
include { BCFTOOLS_RENAME_CHROMS          } from './modules/bcftools/rename_chroms/main'
include { BCFTOOLS_EXTRACT_STRAINS        } from './modules/bcftools/extract_strains/main'
include { BCFTOOLS_CREATE_GENOTYPE_MATRIX } from './modules/bcftools/create_genotype_matrix/main'
include { PLINK_RECODE_VCF as PLINK_RECODE_MS_VCF } from './modules/plink/recode_vcf/main'
include { PLINK_RECODE_VCF as PLINK_RECODE_CV_VCF } from './modules/plink/recode_vcf/main'
include { PLINK_UPDATE_BY_H2              } from './modules/plink/update_by_h2/main'
include { R_FIND_GENOTYPE_MATRIX_EIGEN    } from './modules/r/find_genotype_matrix_eigen/main'
include { PYTHON_SIMULATE_EFFECTS_GLOBAL   } from './modules/python/simulate_effects_global/main'
include { R_SIMULATE_EFFECTS_LOCAL        } from './modules/r/simulate_effects_local/main'
include { R_SIMULATE_EFFECTS_GLOBAL       } from './modules/r/simulate_effects_global/main'
include { R_GET_GCTA_INTERVALS           } from './modules/r/get_gcta_intervals/main'
include { R_ASSESS_SIMS                   } from './modules/r/assess_sims/main'
include { GCTA_SIMULATE_PHENOTYPES        } from './modules/gcta/simulate_phenotypes/main'
// Named aliases for each fanout arm of GCTA_MAKE_GRM. Each alias is a distinct
// process instance in the DAG with its own work dir, logs, and (importantly)
// its own DataflowBroadcast output channels — preventing shared-reference
// issues for further downstream consumers (see ch_grm fanout below).
include { GCTA_MAKE_GRM as GCTA_MAKE_GRM_INBRED } from './modules/gcta/make_grm/main'
include { GCTA_MAKE_GRM as GCTA_MAKE_GRM_LOCO     } from './modules/gcta/make_grm/main'
// Named aliases for the 4-way (mode × type) fanout of GCTA_PERFORM_GWA. Each
// alias is a distinct process instance in the DAG with its own work dir, logs,
// and DataflowBroadcast output channels — preventing shared-FileHolder issues
// at the pca/nopca fanout (same hazard class as the inbred/loco fanout that
// drove the GCTA_MAKE_GRM aliasing above). Outputs are re-mixed at the next
// consumer.
include { GCTA_PERFORM_GWA as GCTA_PERFORM_GWA_INBRED_PCA   } from './modules/gcta/perform_gwa/main'
include { GCTA_PERFORM_GWA as GCTA_PERFORM_GWA_INBRED_NOPCA } from './modules/gcta/perform_gwa/main'
include { GCTA_PERFORM_GWA as GCTA_PERFORM_GWA_LOCO_PCA     } from './modules/gcta/perform_gwa/main'
include { GCTA_PERFORM_GWA as GCTA_PERFORM_GWA_LOCO_NOPCA   } from './modules/gcta/perform_gwa/main'

// Database migration modules
include { DB_MIGRATION_WRITE_MARKER_SET       } from './modules/db_migration/write_marker_set/main'
include { DB_MIGRATION_WRITE_GWA_TO_DB        } from './modules/db_migration/write_gwa_to_db/main'
include { DB_MIGRATION_AGGREGATE_METADATA     } from './modules/db_migration/aggregate_metadata/main'
include { DB_MIGRATION_ANALYZE_QTL            } from './modules/db_migration/analyze_qtl/main'
include { DB_MIGRATION_ASSESS_SIMS            } from './modules/db_migration/assess_sims/main'
include { DB_MIGRATION_WRITE_GENOTYPE_MATRIX  } from './modules/db_migration/write_genotype_matrix/main'
include { DB_MIGRATION_WRITE_TRAIT_DATA       } from './modules/db_migration/write_trait_data/main'

// extractVcfReleaseId: normalizes the strainfile vcf column to a stable 8-digit CaeNDR release date.
// Called from inside a .multiMap{} closure — uses throw, not the NXF error keyword.
// Canonical species allowlist. Single source of truth for both strainfile
// validation and CaeNDR URL construction. Mirrored by .SUPPORTED_SPECIES in
// R/database.R — keep in sync when adding species.
SUPPORTED_SPECIES = ['c_elegans', 'c_briggsae', 'c_tropicalis'] as Set

def extractVcfReleaseId(String vcf) {
    if (vcf ==~ /^\d{8}$/) return vcf
    def matcher = new File(vcf).name =~ /(\d{8})/
    if (matcher.find()) return matcher.group(1)
    throw new IllegalArgumentException(
        "Cannot extract release ID from VCF path '${vcf}' — no 8-digit date found in filename. " +
        "Use YYYYMMDD format in the vcf column or rename the file."
    )
}

// resolveVcf: resolves the strainfile vcf column to an accessible path or URL.
//   - 8-digit date → constructs CaeNDR URL per species
//   - Relative path (does not start with / or http) → prefixed with projectDir
//     (supports test strainfiles that use paths like data/test/test.vcf.gz)
//   - Absolute path or http(s) URL → returned as-is
// Called during channel construction (head-node phase) — errors surface before SLURM submission.
def resolveVcf(String vcf, String species, String projectDir) {
    if (vcf ==~ /^\d{8}$/) {
        if (!SUPPORTED_SPECIES.contains(species)) {
            throw new IllegalArgumentException(
                "Unrecognized species '${species}' for CaeNDR URL construction. " +
                "Supported: ${SUPPORTED_SPECIES}"
            )
        }
        switch (species) {
            case 'c_elegans':    return "https://caendr.org/download/WI.${vcf}.hard-filter.isotype.vcf.gz"
            case 'c_briggsae':   return "https://caendr.org/download/CB.${vcf}.hard-filter.isotype.vcf.gz"
            case 'c_tropicalis': return "https://caendr.org/download/CT.${vcf}.hard-filter.isotype.vcf.gz"
        }
    }
    if (!vcf.startsWith('/') && !vcf.startsWith('http')) {
        return "${projectDir}/${vcf}"  // relative path → absolute (test profile support)
    }
    return vcf  // absolute path passthrough
}

// =====================================================================
// parseManifestTuple — single source of truth for replay manifest types.
//
// The returned 6-tuple MUST match the channel-tuple types produced by
// the pre-process grid construction in main.nf. Mismatched types cause
// Set.contains() to return false for every tuple, silently emptying the
// replay filter and producing a no-op replay run.
//
// maf is intentionally excluded: it is 1:1 with group in the strainfile
// schema, so group alone identifies the marker set. maf still appears as
// a column in replay.tsv and in the channel tuple (it is used in output
// filenames across the pipeline), but it is not part of the Set key.
//
// If a source expression in main.nf changes its emitted type, update this
// helper in the same commit. The runtime assertions in the replay_set
// block will surface drift loudly on the next replay run.
//
// Index map:  0=species  1=group  2=nqtl  3=effect  4=h2  5=rep
// =====================================================================
def parseManifestTuple(Map cols) {
    [
        cols.species,        // String  — matches NF_TRAP_PAYLOAD "species" field
        cols.group,          // String  — matches cleanRow.group
        cols.nqtl,           // String  — matches splitCsv().map { it[0] }
        cols.effect,         // String  — matches splitCsv().map { it[0] }
        cols.h2 as Float,    // Float   — matches h2 channel value after .combine()
        cols.rep as Integer  // Integer — matches Channel.of(1..params.reps)
    ]
}

// Rewrap paths to force fresh FileHolder allocation per arm.
// REQUIRED to avoid race conditions when sibling TaskRuns stage the same
// upstream files concurrently under SLURM array execution. NF caches
// FileHolder instances keyed by Path object identity; sharing a Path
// reference across tuples causes shared FileHolder allocation, which races
// under concurrent staging (NoSuchFileException / truncated stage-ins /
// symlink collisions in work/stage-*). The string round-trip forces
// FileHelper.asPath() to construct a fresh Path per arm.
// Do not remove or simplify without verifying the race no longer occurs at
// scale. Tested with Nextflow 24.10.x (manifest.nextflowVersion guard).
def freshFiles = { paths ->
    if (paths == null) return []
    paths.collect { p ->
        if (p == null) return null
        if (p instanceof List || p instanceof Collection) {
            return p.collect { file(it.toString()) }
        }
        return file(p.toString())
    }
}

workflow {
    main:
    ch_versions = Channel.empty()

    date = new Date().format( 'yyyyMMdd' )

    // Resolve CV (causal variant pool) parameters
    // No 'def' — these become script bindings accessible in workflow.onComplete
    cv_maf = params.cv_maf != null ? params.cv_maf as Float : null
    cv_ld  = params.cv_ld  != null ? params.cv_ld  as Float : 0.8f
    assert cv_ld > 0f && cv_ld < 1.0f : "cv_ld must be in (0, 1), got: ${cv_ld}"
    if (cv_maf != null && cv_maf > 0.5f) {
        error "cv_maf must be in (0, 0.5], got: ${cv_maf} — values above 0.5 are major allele frequency"
    }

    // Set default values for parameters
    if (params.nqtl == null){
        nqtl_file = "${workflow.projectDir}/data/simulate_nqtl.csv"
    } else {
        nqtl_file = params.nqtl
    }
    if (params.h2 == null){
        h2_file = "${workflow.projectDir}/data/simulate_h2.csv"
    } else {
        h2_file = params.h2
    }
    if (params.effect == null){
        effect_file = "${workflow.projectDir}/data/simulate_effect_sizes.csv"
    } else {
        effect_file = params.effect
    }
    if (params.strainfile == null){
        strainfile = "${workflow.projectDir}/data/test_strains.txt"
    } else {
        strainfile = file(params.strainfile).toAbsolutePath().toString()
    }
    ch_input_files = Channel.of(
        file(strainfile),
        file(nqtl_file),
        file(h2_file),
        file(effect_file)
    )
    if (params.mito_name == null){
        mito_name = "MtDNA"
    } else {
        mito_name = params.mito_name
    }
    if (params.simulate_qtlloc == null){
        simulate_qtlloc = false
    } else {
        simulate_qtlloc = params.simulate_qtlloc
    }

    // set help message
    if (params.help) {
        log.info '''
        '''
        log.info "----------------------------------------------------------------"
        log.info "                      USAGE                                     "
        log.info "----------------------------------------------------------------"
        log.info " "
        log.info "nextflow andersenlab/nemascan-sim-nf --strainfile /path/to/strainfile --vcf /path/to/vcf -output-dir my-results"
        log.info " "
        log.info "Profiles available:"
        log.info "rockfish              Profile        Perform selected analysis on Rockfish (default GWA mapping)"
        log.info " "
        log.info "Mandatory argument (General):"
        log.info "--strainfile      File               A tab-separated file with header: group, species, vcf, ms_maf, ms_ld, strains. One row per strain group. vcf accepts a CaeNDR release date (8-digit) or an absolute path."
        log.info " "
        log.info "Optional arguments (General):"
        log.info "--nqtl            File               A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: data/simulate_nqtl.csv)"
        log.info "--h2              File               A CSV file with phenotype heritability, one value per line (Default is located: data/simulate_h2.csv)"
        log.info "--reps             Integer            The number of replicates to simulate per number of QTL and heritability (Default: 2)"
        log.info "--cv_maf          Decimal            Minor allele frequency threshold for causal variant pool (Default: per-group ms_maf from strainfile)"
        log.info "--cv_ld           Decimal            LD R² threshold for causal variant pool pruning (Default: 0.8)"
        log.info "--effect          File               A CSV file where each line is an effect size range (e.g. 0.2-0.3) to test for simulations (Default: data/simulate_effect_sizes.csv)"
        log.info "--qtlloc          File               A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])"
        log.info "--sthresh         String             Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
        log.info "--alpha           Decimal            Significance level for Bonferroni and EIGEN threshold calculation (Default: 0.05)"
        log.info "--group_qtl       Integer            If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
        log.info "--ci_size         Integer            Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
        log.info "--sparse_cut      Decimal            Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
        log.info "--mito_name       Strain             Name of mitochondrial chromosome"
        log.info "--simulate_qtlloc Boolean            Whether to simulate QTLs in specific genomic regions (Default: false)"
        log.info "--legacy_assess   Boolean            Run legacy R-based QTL detection in parallel with DB path for cross-validation (Default: false)"
        log.info "--output_dir      String             Output directory name (Default: Analysis_Results-{date}). Also settable via Nextflow native -output-dir flag."
        log.info " "


        exit 1
    } else { // set log info
        log.info '''


    '''
        log.info ""
        log.info "Strainfile                              = ${strainfile}"
        log.info "Number of QTLs/phenotype simulated      = ${nqtl_file}"
        log.info "Phenotype heritability file             = ${h2_file}"
        log.info "Number of replicates to simulate        = ${params.reps}"
        log.info "Causal variant MAF                      = ${cv_maf ?: '(per-group ms_maf)'}"
        log.info "Causal variant LD threshold             = ${cv_ld}"
        log.info "Effect size range file                  = ${effect_file}"
        log.info "Genome range file                       = ${params.qtlloc}"
        log.info "Significance Thresholds                 = BF, EIGEN"
        log.info "Window for combining QTLs               = ${params.group_qtl}"
        log.info "Number of SNVs to define QTL CI         = ${params.ci_size}"
        log.info "Relatedness cutoff                      = ${params.sparse_cut}"
        log.info "Mitochondrial chromosome name           = ${mito_name}"
        log.info "Simulate QTLs in specific regions       = ${simulate_qtlloc}"
        log.info "Output directory                        = ${workflow.outputDir}"
        log.info ""
    }

    // Pre-create the failures directory before any task can launch.
    // workflow.onStart does not exist in NF v24 — use the workflow body directly.
    def failuresDir = file("${workflow.outputDir}/.failures")
    failuresDir.mkdirs()
    assert failuresDir.exists() :
        "Failed to create ${failuresDir} — check permissions on the output directory"

    // Replay mode: parse the manifest once at workflow start, build a typed Set.
    // null when --replay is not set (normal runs — no filtering applied).
    def replay_set = params.replay
        ? {
              def lines = file(params.replay).readLines()
              assert lines && lines[0].startsWith('session\t') :
                  "replay.tsv header does not start with 'session\\t' — file may be malformed or missing its header"
              def header = lines[0].split('\t')
              def required = ['species', 'group', 'nqtl', 'effect', 'h2', 'rep'] as Set
              assert required.every { header.contains(it) } :
                  "replay.tsv header is missing required columns; found: ${header.toList()} — need: ${required}"
              def tuples = lines.drop(1).collect { line ->
                  def fields = line.split('\t')
                  assert fields.size() == header.size() :
                      "row tokenization failure — column count mismatch: ${line}"
                  def row = [header, fields].transpose().collectEntries()
                  parseManifestTuple(row)
              }.toSet()

              // Cheap type-alignment assertions — surface coercion drift before the filter runs.
              // Tuple index map: 0=species 1=group 2=nqtl 3=effect 4=h2 5=rep
              if (tuples) {
                  def s = tuples.first()
                  assert s[0].class == String  : "replay_set species type is ${s[0].class} — expected String; check parseManifestTuple vs main.nf param-grid"
                  assert s[1].class == String  : "replay_set group type is ${s[1].class} — expected String; check parseManifestTuple vs main.nf param-grid"
                  assert s[2].class == String  : "replay_set nqtl type is ${s[2].class} — expected String; check parseManifestTuple vs main.nf param-grid"
                  assert s[4].class == Float   : "replay_set h2 type is ${s[4].class} — expected Float; check parseManifestTuple vs main.nf param-grid"
                  assert s[5].class == Integer : "replay_set rep type is ${s[5].class} — expected Integer; check parseManifestTuple vs main.nf param-grid"
              }

              log.info "Replay mode: filtering channels to ${tuples.size()} slot(s) from ${params.replay}"
              tuples
          }()
        : null

    // Parse 6-column strainfile and fan out per-row parameters.
    //
    // Species is normalized (lowercased + trimmed) and validated against
    // SUPPORTED_SPECIES once, here. From this point on the canonical species
    // value lives in `meta.species` and is propagated through every per-group
    // sub-channel. No downstream code should re-read or re-normalize species
    // from the strainfile row.
    ch_strain_sets = Channel.fromPath(strainfile)
        .splitCsv(sep: "\t", header: true)
        .map { row ->
            // M8 fix: normalize whitespace to handle Windows CRLF line endings
            def cleanRow = row.collectEntries { k, v -> [k.trim(), v?.trim()] }
            def speciesNorm = cleanRow.species?.toLowerCase()
            if (!SUPPORTED_SPECIES.contains(speciesNorm)) {
                error "Unknown species '${cleanRow.species}' in strainfile row for group '${cleanRow.group}' — supported: ${SUPPORTED_SPECIES}"
            }
            if (cleanRow.ms_maf.toFloat() <= 0 || cleanRow.ms_maf.toFloat() > 0.5) {
                error "ms_maf '${cleanRow.ms_maf}' out of range (0, 0.5] in group '${cleanRow.group}' — values above 0.5 are major allele frequency, not MAF"
            }
            if (cleanRow.ms_ld.toFloat() <= 0 || cleanRow.ms_ld.toFloat() >= 1.0) {
                error "ms_ld '${cleanRow.ms_ld}' out of range (0, 1) in group '${cleanRow.group}'"
            }
            [
                [id: cleanRow.group, species: speciesNorm],
                cleanRow.vcf,
                cleanRow.ms_maf.toFloat(),
                cleanRow.ms_ld.toFloat(),
                cleanRow.strains
            ]
        }

    ch_strain_sets
        .multiMap { meta, vcf, ms_maf, ms_ld, strains ->
            def cv_maf_eff = cv_maf != null ? cv_maf : ms_maf
            if (cv_maf != null && cv_maf > ms_maf) {
                log.warn "cv_maf (${cv_maf}) > ms_maf (${ms_maf}) for group ${meta.id}: " +
                         "the CV pool is a strict subset of the marker SNP set. " +
                         "This is likely unintentional."
            }
            marker_set_params:  [meta.id, ms_maf, meta.species, extractVcfReleaseId(vcf), ms_ld, strains, strainfile]
            vcf_per_group:      [meta, vcf, strains]
            ms_maf_vals:        [meta, ms_maf]
            ms_ld_vals:         [meta, ms_ld]
            cv_maf_vals:        [meta, cv_maf_eff]
            species_per_group:  [meta.id, meta.species]
        }
        .set { ch_sf }

    // Single attachment point for species: combined onto the per-cell sim grid
    // exactly once. From PYTHON_SIMULATE_EFFECTS_GLOBAL onward, species rides
    // along inside each process's params output tuple, so downstream consumers
    // get species "for free" without additional forks or combines.
    //
    // ch_sf.species_per_group is consumed exactly once, here. If you ever need
    // another standalone copy elsewhere, restore a .tap{} fork — do not consume
    // the queue channel twice directly.
    ch_sf.species_per_group
        .set { ch_species_for_grid }

    // G1: Fork cv_maf_vals for four consumers: PLINK_RECODE_CV_VCF, write_trait_data,
    // analysis params, and write_gwa_to_db. A single multiMap sub-channel cannot be
    // consumed by more than one operator; .tap{} creates a broadcast fork.
    ch_sf.cv_maf_vals
        .tap { ch_cv_maf_for_plink }        // → PLINK_RECODE_CV_VCF
        .tap { ch_cv_maf_for_trait }        // → combine in merge chain (G3)
        .tap { ch_cv_maf_for_gwa_write }    // → extend ch_gwa_db_inputs (G5)
        .map { meta, cv_maf_eff -> tuple(meta.id, cv_maf_eff) }
        .set { ch_cv_maf_keyed_for_analysis }  // → extend analysis params (G4)

    ch_cv_maf_keyed_for_trait = ch_cv_maf_for_trait
        .map { meta, cv_maf_eff -> tuple(meta.id, cv_maf_eff) }

    ch_cv_maf_keyed_for_gwa_write = ch_cv_maf_for_gwa_write
        .map { meta, cv_maf_eff -> tuple(meta.id, cv_maf_eff) }

    // ch_sf.marker_set_params is a queue channel — fan out with .tap{} before
    // 2 subscribers consume it (write_marker_set, write_genotype_matrix).
    // Without this, DSL2 distributes emissions round-robin among competing readers.
    // write_gwa_to_db no longer needs a fork here; it reads species/vcf_release_id/ms_ld
    // from marker_set_metadata.parquet at runtime via read_marker_set_metadata().
    ch_sf.marker_set_params
        .tap { ch_marker_set_params_for_ms }
        .set { ch_marker_set_params_for_gm }

    ch_vcf_per_group = ch_sf.vcf_per_group
        .map { meta, vcf, strains ->
            def vcf_path = resolveVcf(vcf, meta.species, workflow.projectDir.toString())
            [meta, file(vcf_path), file("${vcf_path}.tbi"), strains]
        }

    // Fan out ch_vcf_per_group to separate sub-channels (avoids queue channel consumption)
    ch_vcf_per_group
        .multiMap { meta, vcf, tbi, strains ->
            for_contig:  [meta, vcf, tbi]
            for_extract: [meta, vcf, tbi]
            for_strains: [meta, strains]
        }
        .set { ch_vcf_groups }

    // Get contig data from VCF file
    // .first() — one representative VCF suffices; all supported species share same chromosome names
    LOCAL_GET_CONTIG_INFO( ch_vcf_groups.for_contig.first() )
    ch_mito_num = LOCAL_GET_CONTIG_INFO.out.mapping.splitCsv(sep:"\t")
        .filter{ row -> row[0] == mito_name }
        .map{ row -> row[1] }
        .first()

    BCFTOOLS_EXTRACT_STRAINS( ch_vcf_groups.for_extract, ch_vcf_groups.for_strains )
    ch_versions = ch_versions.mix(BCFTOOLS_EXTRACT_STRAINS.out.versions)

    // Extract desired strain sets
    BCFTOOLS_RENAME_CHROMS(
        BCFTOOLS_EXTRACT_STRAINS.out.vcf,
        LOCAL_GET_CONTIG_INFO.out.mapping
        )
    ch_versions = ch_versions.mix(BCFTOOLS_RENAME_CHROMS.out.versions)

    // Recode the VCF file and create plink formatted files
    // Fork renamed VCF before two parallel PLINK calls (avoids queue channel consumption)
    BCFTOOLS_RENAME_CHROMS.out.vcf
        .tap { ch_renamed_for_ms }
        .tap { ch_renamed_for_cv }

    // Marker SNP selection — uses per-group ms_maf and ms_ld from strainfile
    PLINK_RECODE_MS_VCF(
        ch_renamed_for_ms,
        ch_mito_num,
        ch_sf.ms_maf_vals.map { meta, ms_maf -> ms_maf },
        ch_sf.ms_ld_vals.map  { meta, ms_ld  -> ms_ld  }
        )
    ch_versions = ch_versions.mix(PLINK_RECODE_MS_VCF.out.versions)

    // Causal variant pool selection — uses global cv_maf/cv_ld (or per-group ms_maf fallback)
    PLINK_RECODE_CV_VCF(
        ch_renamed_for_cv,
        ch_mito_num,
        ch_cv_maf_for_plink.map { meta, cv_maf_eff -> cv_maf_eff },
        Channel.value(cv_ld)
        )
    ch_versions = ch_versions.mix(PLINK_RECODE_CV_VCF.out.versions)

    // Create plaintext genotype matrix
    BCFTOOLS_CREATE_GENOTYPE_MATRIX(
        PLINK_RECODE_MS_VCF.out.vcf,
        PLINK_RECODE_MS_VCF.out.markers
        )
    ch_versions = ch_versions.mix(BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.versions)

    // Find eigen values for genotype matrix
    ch_chrom_nums = LOCAL_GET_CONTIG_INFO.out.mapping
        .splitCsv(sep: "\t")
        .filter { it -> it[0] != mito_name }
        .map { it -> it[1] }
        .toSortedList()

    R_FIND_GENOTYPE_MATRIX_EIGEN(
        BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix,
        Channel.fromPath("${workflow.projectDir}/bin/Get_GenoMatrix_Eigen.R").first(),
        ch_chrom_nums
        )
    ch_versions = ch_versions.mix(R_FIND_GENOTYPE_MATRIX_EIGEN.out.versions)

    // Concatenate eigen files with plink and genotype matrix sets
    ch_eigens = R_FIND_GENOTYPE_MATRIX_EIGEN.out.eigen
        .groupTuple(by: [0, 1])
        .map{ it: [it[0], it[1], it[3]] }

    LOCAL_COMPILE_EIGENS( ch_eigens )

    // CV binary — drop maf (cv_maf_eff) so combine(by:0) keys on group only
    ch_cv_plink_for_combine = PLINK_RECODE_CV_VCF.out.plink
        .map { group, _maf, bed, bim, fam, map_f, nosex, ped, log_f ->
            [group, bed, bim, fam, map_f, nosex, ped, log_f]
        }

    // Compile required files for simulations by strain group and MAF
    // Extends with CV binary so PYTHON_SIMULATE_EFFECTS_GLOBAL can sample non-marker causal variants
    ch_plink_genomat_eigen = PLINK_RECODE_MS_VCF.out.plink
        .join(BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix, by: [0, 1])
        .join(LOCAL_COMPILE_EIGENS.out.tests, by: [0, 1])
        .combine(ch_cv_plink_for_combine, by: 0)
    // Resulting tuple (18 elements):
    //   [group, ms_maf, ms_bed, ms_bim, ms_fam, ms_map, ms_nosex, ms_ped, ms_log,
    //    gm, n_indep_tests,
    //    cv_bed, cv_bim, cv_fam, cv_map, cv_nosex, cv_ped, cv_log]
    
    // // Simulate QTL or genome
    // if (simulate_qtlloc){
    //     R_SIMULATE_EFFECTS_LOCAL( ch_plink_genomat_eigen,
    //                               Channel.fromPath(params.qtlloc),
    //                               Channel.fromPath("${workflow.projectDir}/bin/Create_Causal_QTLs.R").first(),
    //                               Channel.of(1..params.reps).toSortedList(),
    //                               Channel.fromPath(nqtl_file)
    //                                 .splitCsv()
    //                                 .map{ it: it[0] }
    //                                 .toSortedList(),
    //                               Channel.fromPath(effect_file)
    //                                 .splitCsv()
    //                                 .map{ it: it[0] }
    //                                 .toSortedList() 
    //                             )
    //     ch_versions = ch_versions.mix(R_SIMULATE_EFFECTS_LOCAL.out.versions)
    //     ch_sim_phenos = R_SIMULATE_EFFECTS_LOCAL.out.causal
    //     ch_sim_plink = R_SIMULATE_EFFECTS_LOCAL.out.plink
    // } else {

    def replay_predicate = { args -> replay_set.contains([args[0], args[1], args[3], args[4], args[5], args[6]]) }

    // Build the full (species, group, maf, nqtl, effect, h2, rep, ...files...) grid
    // at channel level before PYTHON_SIMULATE_EFFECTS_GLOBAL.
    // Replaces: each rep (inside process), each nqtl (inside process), each effect (inside process).
    // GCTA_SIMULATE_PHENOTYPES will also lose each h2 — h2 flows through via the causal output.
    // ch_sim_plink / ch_sim_cv_plink are NOT filtered — keyed only on (group, maf).
    ch_sim_grid = ch_plink_genomat_eigen
        .combine(ch_species_for_grid, by: 0)
        // → (group, maf, ms_bed…ms_log, gm, n_indep_tests, cv_bed…cv_log, species)
        .combine(Channel.fromPath(nqtl_file).splitCsv().map { it[0] })
        .combine(Channel.fromPath(effect_file).splitCsv().map { it[0] })
        .combine(Channel.fromPath(h2_file).splitCsv().map { it[0] })
        .combine(Channel.of(1..params.reps))
        // → reorder to canonical (species, group, maf, nqtl, effect, h2, rep, ...files...)
        .map { group, maf,
               ms_bed, ms_bim, ms_fam, ms_map, ms_nosex, ms_ped, ms_log,
               gm, n_indep_tests,
               cv_bed, cv_bim, cv_fam, cv_map, cv_nosex, cv_ped, cv_log,
               species, nqtl, effect, h2, rep ->
            tuple(species, group, maf, nqtl, effect, h2 as Float, rep,
                  ms_bed, ms_bim, ms_fam, ms_map, ms_nosex, ms_ped, ms_log,
                  gm, n_indep_tests,
                  cv_bed, cv_bim, cv_fam, cv_map, cv_nosex, cv_ped, cv_log)
        }
        .filter { params.replay ? replay_predicate(it) : true }


    PYTHON_SIMULATE_EFFECTS_GLOBAL(
        ch_sim_grid,
        Channel.fromPath("${workflow.projectDir}/bin/create_causal_vars.py").first()
        )
    ch_versions = ch_versions.mix(PYTHON_SIMULATE_EFFECTS_GLOBAL.out.versions)
    ch_sim_phenos    = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.causal
    ch_sim_plink     = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.plink
    ch_sim_cv_plink  = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.cv_plink
    // }

    // G2: causal_genotypes emits (group, maf, nqtl, effect, rep, h2, geno_file) — h2 is already
    // in the tuple (fanned out at channel level before the process). No .combine(h2) needed.
    // The join by [0,1,2,3,4,5] in the merge chain (G3) keys on (group, maf, nqtl, effect, rep, h2).
    ch_causal_geno_fanned = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.causal_genotypes

    // Cleanup gate — deletes failed-slot DB files before any trait-simulation task launches.
    // DB_CLEAN_REPLAY_SLOTS runs once per --replay invocation; its done sentinel is a value
    // channel (.first()) so it broadcasts to every downstream tuple rather than consuming once.
    // ch_sim_plink / ch_sim_cv_plink are gated here even though they are not filtered (their
    // keys are group/maf only), to prevent any GWA write from starting before cleanup finishes.
    def db_output_dir = params.db_output ? file(params.db_output).toAbsolutePath().toString() : file("${workflow.outputDir}/db").toAbsolutePath().toString()
    log.info "Database output directory: ${db_output_dir}"

    def ch_cleanup_done
    if (params.replay) {
        DB_CLEAN_REPLAY_SLOTS(file(params.replay), db_output_dir, cv_maf, cv_ld)
        ch_cleanup_done = DB_CLEAN_REPLAY_SLOTS.out.done.first()
    } else {
        ch_cleanup_done = Channel.value('no_cleanup_needed')
    }

    ch_sim_phenos   = ch_sim_phenos.combine(ch_cleanup_done).map   { it[0..-2] }
    ch_sim_plink    = ch_sim_plink.combine(ch_cleanup_done).map    { it[0..-2] }
    ch_sim_cv_plink = ch_sim_cv_plink.combine(ch_cleanup_done).map { it[0..-2] }

    // PYTHON_SIMULATE_EFFECTS_GLOBAL emits species at the end of its .causal
    // tuple (see modules/python/simulate_effects_global/main.nf). ch_sim_phenos
    // shape: (group, maf, nqtl, effect, rep, h2, causal_file, species) —
    // matches GCTA_SIMULATE_PHENOTYPES input directly. No combine needed.
    //
    // Simulate phenotypes using the CV binary (--bfile CV_TO_SIMS) so GCTA can locate
    // non-marker causal variant genotypes. The MS binary (ch_sim_plink) passes through
    // to downstream GWA mapping unchanged.
    GCTA_SIMULATE_PHENOTYPES(
        ch_sim_phenos,
        ch_sim_cv_plink,   // CV binary — used as --bfile for GCTA simulation
        ch_sim_plink       // MS binary — passed through to downstream GWA
        )
    ch_versions = ch_versions.mix(GCTA_SIMULATE_PHENOTYPES.out.versions)

    // Fork GCTA_SIMULATE_PHENOTYPES outputs for two consumers:
    //   1. .merge() chain → ch_trait (DB trait data writes)
    //   2. PLINK_UPDATE_BY_H2 (downstream GWA pipeline)
    // Without explicit .tap{}, implicit broadcast shares mutable ArrayList
    // references between consumers. Under high concurrency with SLURM job
    // arrays + Singularity, BashWrapperBuilder.createContainerBuilder can
    // race with the merge chain thread -> ConcurrentModificationException.
    GCTA_SIMULATE_PHENOTYPES.out.params
        .tap { ch_gcta_params_for_trait }
        .set { ch_gcta_params_for_plink }

    GCTA_SIMULATE_PHENOTYPES.out.pheno
        .tap { ch_gcta_pheno_for_trait }
        .set { ch_gcta_pheno_for_plink }

    GCTA_SIMULATE_PHENOTYPES.out.plink
        .tap { ch_gcta_plink_for_trait }
        .set { ch_gcta_plink_for_plink }

    // -- TRAIT DATA WRITES -------------------------------------------------------
    // Source: GCTA_SIMULATE_PHENOTYPES (pre-upscaled phenotype, pre-mode-crossing)
    //
    // GCTA_SIMULATE_PHENOTYPES runs 1x per trait (before mode crossing at
    // GCTA_MAKE_GRM line 260). No mode deduplication needed -- each emission
    // is a unique trait.
    //
    // Trait data (metadata, causal variants, phenotype, and causal genotypes) are
    // stored in the DB for future analysis. These data are NOT consumed by
    // ASSESS_SIMS -- no barrier is needed.

    // G3: Combine with cv_maf and join with causal genotypes before multiMap so that
    // write_pheno/write_par stay emission-order-aligned with write_params.
    //
    // After .combine(ch_cv_maf_keyed_for_trait, by: 0): 21-element tuple (cv_maf_eff at [20])
    // After .join(ch_causal_geno_fanned, by: [0..5]):   22-element tuple (geno_file at [21])
    //
    // WARNING: .merge() aligns channels by emission order, NOT by key matching.
    // Do NOT insert .filter(), .map(), .branch(), or any reordering operator
    // between GCTA_SIMULATE_PHENOTYPES.out.* and this .merge() chain — doing
    // so will silently misalign phenotype/causal files with metadata.
    // Runtime validation in write_trait_data.R catches misalignment.
    // .combine(by:0) and .join(by:[0..5]) are key-based (order-preserving) and
    // therefore safe per the CLAUDE.md merge ordering constraint.
    //
    // Merged tuple structure (22 elements):
    //   [0-5]   params:   group, maf, nqtl, effect, rep, h2
    //   [6]     params:   species  (forwarded from GCTA_SIMULATE_PHENOTYPES.out.params)
    //   [7-8]   pheno:    phen_path, par_path
    //   [9-19]  plink:    group, maf, bed, bim, fam, map, nosex, ped, log, gm, n_indep_tests
    //   [20]    cv_maf:   cv_maf_effective (from combine)
    //   [21]    geno:     causal_genotypes TSV path (from join)
    //
    // write_trait_data.R reads species from marker_set_metadata.parquet at
    // runtime; the species value in this tuple is not consumed by the trait
    // write — it's forwarded only so downstream processes that DO need it
    // (PLINK_UPDATE_BY_H2 → GCTA_MAKE_GRM → GCTA_PERFORM_GWA) can read it
    // from ch_gcta_params_for_plink without a fresh combine.
    //
    // IMPORTANT: If GCTA_SIMULATE_PHENOTYPES output channels change, update
    // these indices.
    ch_gcta_params_for_trait
        .merge(ch_gcta_pheno_for_trait)
        .merge(ch_gcta_plink_for_trait)
        .combine(ch_cv_maf_keyed_for_trait, by: 0)
        .join(ch_causal_geno_fanned, by: [0, 1, 2, 3, 4, 5])
        .multiMap { it ->
            assert it.size() == 22 : "Expected 22-element tuple, got ${it.size()}. First elements: ${it.take(7)}. Check output channels in modules/gcta/simulate_phenotypes/main.nf"
            // trait_id is computed in R by write_trait_data.R via generate_trait_id()
            // — no cross-language hash computation (addresses review B2)
            write_params:  tuple(it[0], it[1], it[2], it[3], it[4], it[5])
            write_pheno:   it[7]   // pre-upscaled .phen file
            write_par:     it[8]   // causal variant .par file
            cv_maf:        it[20]  // cv_maf_effective (from combine with ch_cv_maf_keyed_for_trait)
            causal_geno:   it[21]  // causal genotype TSV (from join with ch_causal_geno_fanned)
        }
        .set { ch_trait }

    // WRITE_TRAIT_DATA is gated below (after ch_marker_barrier) because
    // write_trait_data.R reads marker_set_metadata.parquet to resolve the
    // marker set ID. Without the gate it races WRITE_MARKER_SET.

    // Update plink data by heritability
    PLINK_UPDATE_BY_H2(
        ch_gcta_params_for_plink,
        ch_gcta_plink_for_plink,
        ch_gcta_pheno_for_plink
        )
    ch_versions = ch_versions.mix(PLINK_UPDATE_BY_H2.out.versions)

    // Create genetic relatedness matrix.
    // PLINK_UPDATE_BY_H2 emits a single combined tuple (params + plink + pheno)
    // per task. We fan out across mode (inbred/loco) via flatMap, rewrapping each
    // path with freshFiles() to force NF to allocate independent FileHolder
    // instances per arm. Each arm is then routed to its own named alias of
    // GCTA_MAKE_GRM (GCTA_MAKE_GRM_INBRED / GCTA_MAKE_GRM_LOCO) so
    // sibling TaskRuns land in distinct work-dir subtrees with distinguishable
    // process names in `nextflow log` and traces.
    //
    // freshFiles() rewrapping is REQUIRED — forces fresh FileHolder per arm to
    // avoid SLURM-array staging race (NoSuchFileException, truncated stage-ins,
    // symlink collisions). Do not refactor back to .combine() or remove the
    // file(p.toString()) round-trip without verifying the race at scale.
    // See issues/175-gcta-concurrent-mod-exception/.
    // PLINK_UPDATE_BY_H2.out.out tuple shape (18 elements):
    //   group, maf, nqtl, effect, rep, h2,
    //   bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
    //   pheno, par,
    //   species  (appended in modules/plink/update_by_h2/main.nf)
    //
    // Species rides at the end of the upstream tuple; we destructure it and
    // include it in the params slot for GCTA_MAKE_GRM (which expects species
    // at the end of its params input).
    PLINK_UPDATE_BY_H2.out.out
        .flatMap { group, maf, nqtl, effect, rep, h2,
                   bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                   pheno, par, species ->
            def plink_files = [bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests]
            def pheno_files = [pheno, par]
            def prefix      = [group, maf, nqtl, effect, rep, h2]
            [
                // Call freshFiles() SEPARATELY for each arm — sharing a single
                // freshFiles() result list across arms reintroduces the race.
                prefix + ["inbred", "fastGWA", species] + freshFiles(plink_files) + freshFiles(pheno_files),
                prefix + ["loco",   "mlma",    species] + freshFiles(plink_files) + freshFiles(pheno_files)
            ]
        }
        .branch { group, maf, nqtl, effect, rep, h2, mode, suffix, species,
                  bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                  pheno, par ->
            inbred: mode == "inbred"
            loco:   mode == "loco"
        }
        .set { ch_grm_routed }

    ch_grm_routed.inbred
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, species,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par ->
            params:    tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, species)
            plink:     tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno_in:  tuple(pheno, par)
        }
        .set { ch_grm_inbred }

    ch_grm_routed.loco
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, species,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par ->
            params:    tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, species)
            plink:     tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno_in:  tuple(pheno, par)
        }
        .set { ch_grm_loco }

    GCTA_MAKE_GRM_INBRED(
        ch_grm_inbred.params,
        ch_grm_inbred.plink,
        ch_grm_inbred.pheno_in
        )
    GCTA_MAKE_GRM_LOCO(
        ch_grm_loco.params,
        ch_grm_loco.plink,
        ch_grm_loco.pheno_in
        )
    ch_versions = ch_versions
        .mix(GCTA_MAKE_GRM_INBRED.out.versions)
        .mix(GCTA_MAKE_GRM_LOCO.out.versions)

    // Simulate GWA using output from GCTA_MAKE_GRM.
    //
    // Strategy: each alias emits a single combined output tuple, so the two
    // GRM-alias streams are mixed directly — no per-alias reassembly needed.
    // Species rides inside GCTA_MAKE_GRM.out.combined at position [8], so no
    // separate combine is needed — it flows forward into GCTA_PERFORM_GWA.
    //
    // A flatMap with freshFiles fans out pca/nopca arms — same race-defense
    // rewrap as the inbred/loco fanout. The branch routes each (mode, type)
    // tuple to its own multiMap+alias invocation, producing four distinct
    // GCTA_PERFORM_GWA TaskRun streams.
    //
    // freshFiles() is REQUIRED here for the pca/nopca arms — both arms
    // would otherwise share FileHolder references to the same GRM-alias
    // output paths (grm + plink + pheno), reintroducing the SLURM-array
    // staging race that motivated this aliasing.
    //
    // GCTA_MAKE_GRM.out.combined shape (23 elements):
    //   [0-8]   vals:  group, maf, nqtl, effect, rep, h2, mode, suffix, species
    //   [9-11]  grm:   grm.bin, grm.N.bin, grm.id
    //   [12-20] plink: bed, bim, fam, map, nosex, ped, log, gm, n_indep_tests
    //   [21-22] pheno: pheno, par
    ch_grm_inbred_out = GCTA_MAKE_GRM_INBRED.out.combined
    ch_grm_loco_out   = GCTA_MAKE_GRM_LOCO.out.combined

    ch_grm_inbred_out
        .mix(ch_grm_loco_out)
        .flatMap { group, maf, nqtl, effect, rep, h2, mode, suffix, species,
                   grm_bin, grm_n, grm_id,
                   bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                   pheno, par ->
            def prefix      = [group, maf, nqtl, effect, rep, h2, mode, suffix]
            def grm_files   = [grm_bin, grm_n, grm_id]
            def plink_files = [bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests]
            def pheno_files = [pheno, par]
            [
                // freshFiles() called SEPARATELY per arm — sharing a single
                // result list across arms reintroduces the race.
                // Species emitted twice (once per pca/nopca arm) — each arm
                // is an independent task with its own input tuple.
                prefix + ["pca", species]   + freshFiles(grm_files) + freshFiles(plink_files) + freshFiles(pheno_files),
                prefix + ["nopca", species] + freshFiles(grm_files) + freshFiles(plink_files) + freshFiles(pheno_files)
            ]
        }
        .branch { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                  grm_bin, grm_n, grm_id,
                  bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                  pheno, par ->
            inbred_pca:   mode == "inbred" && type == "pca"
            inbred_nopca: mode == "inbred" && type == "nopca"
            loco_pca:     mode == "loco"   && type == "pca"
            loco_nopca:   mode == "loco"   && type == "nopca"
        }
        .set { ch_gwa_routed }

    // Per-arm multiMap — sub-channel `pheno_in` to disambiguate from upstream
    // `.pheno` channel name. Each multiMap constructs fresh tuple() values so
    // sibling alias TaskRuns do not share outer ArrayList references.
    ch_gwa_routed.inbred_pca
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                    grm_bin, grm_n, grm_id,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par ->
            params:    tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species)
            grm:       tuple(grm_bin, grm_n, grm_id)
            plink:     tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno_in:  tuple(pheno, par)
        }
        .set { ch_gwa_inbred_pca }

    ch_gwa_routed.inbred_nopca
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                    grm_bin, grm_n, grm_id,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par ->
            params:    tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species)
            grm:       tuple(grm_bin, grm_n, grm_id)
            plink:     tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno_in:  tuple(pheno, par)
        }
        .set { ch_gwa_inbred_nopca }

    ch_gwa_routed.loco_pca
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                    grm_bin, grm_n, grm_id,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par ->
            params:    tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species)
            grm:       tuple(grm_bin, grm_n, grm_id)
            plink:     tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno_in:  tuple(pheno, par)
        }
        .set { ch_gwa_loco_pca }

    ch_gwa_routed.loco_nopca
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                    grm_bin, grm_n, grm_id,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par ->
            params:    tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species)
            grm:       tuple(grm_bin, grm_n, grm_id)
            plink:     tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno_in:  tuple(pheno, par)
        }
        .set { ch_gwa_loco_nopca }

    GCTA_PERFORM_GWA_INBRED_PCA(
        ch_gwa_inbred_pca.params,
        ch_gwa_inbred_pca.grm,
        ch_gwa_inbred_pca.plink,
        ch_gwa_inbred_pca.pheno_in,
        params.sparse_cut
        )
    GCTA_PERFORM_GWA_INBRED_NOPCA(
        ch_gwa_inbred_nopca.params,
        ch_gwa_inbred_nopca.grm,
        ch_gwa_inbred_nopca.plink,
        ch_gwa_inbred_nopca.pheno_in,
        params.sparse_cut
        )
    GCTA_PERFORM_GWA_LOCO_PCA(
        ch_gwa_loco_pca.params,
        ch_gwa_loco_pca.grm,
        ch_gwa_loco_pca.plink,
        ch_gwa_loco_pca.pheno_in,
        params.sparse_cut
        )
    GCTA_PERFORM_GWA_LOCO_NOPCA(
        ch_gwa_loco_nopca.params,
        ch_gwa_loco_nopca.grm,
        ch_gwa_loco_nopca.plink,
        ch_gwa_loco_nopca.pheno_in,
        params.sparse_cut
        )
    ch_versions = ch_versions
        .mix(GCTA_PERFORM_GWA_INBRED_PCA.out.versions)
        .mix(GCTA_PERFORM_GWA_INBRED_NOPCA.out.versions)
        .mix(GCTA_PERFORM_GWA_LOCO_PCA.out.versions)
        .mix(GCTA_PERFORM_GWA_LOCO_NOPCA.out.versions)

    // ── RE-MIX AT THE NEXT CONSUMER ──────────────────────────────────────
    // Each alias emits a single combined 25-element output tuple. Mix across
    // the 4 aliases collapses the streams, then multiMap breaks out the named
    // sub-channels (params, grm, plink, pheno, gwa) that downstream consumers
    // expect. Because multiMap emits all sub-channels from the same source
    // tuple in lockstep, they remain correctly aligned after the alias split.
    //
    // GCTA_PERFORM_GWA.out.combined shape (25 elements):
    //   [0-9]   vals:  group, maf, nqtl, effect, rep, h2, mode, suffix, type, species
    //   [10-12] grm:   grm.bin, grm.N.bin, grm.id
    //   [13-21] plink: bed, bim, fam, map, nosex, ped, log, gm, n_indep_tests
    //   [22-23] pheno: pheno, par
    //   [24]    gwa:   mapping result file
    ch_inbred_pca_out   = GCTA_PERFORM_GWA_INBRED_PCA.out.combined
    ch_inbred_nopca_out = GCTA_PERFORM_GWA_INBRED_NOPCA.out.combined
    ch_loco_pca_out     = GCTA_PERFORM_GWA_LOCO_PCA.out.combined
    ch_loco_nopca_out   = GCTA_PERFORM_GWA_LOCO_NOPCA.out.combined

    ch_inbred_pca_out
        .mix(ch_inbred_nopca_out)
        .mix(ch_loco_pca_out)
        .mix(ch_loco_nopca_out)
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                    grm_bin, grm_n, grm_id,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par, gwa ->
            params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species)
            grm:    tuple(grm_bin, grm_n, grm_id)
            plink:  tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno:  tuple(pheno, par)
            gwa:    gwa
        }
        .set { ch_perform_gwa }

    // Fork GCTA_PERFORM_GWA outputs for multiple consumers:
    //   .params → DB write path, DB analysis path, [legacy intervals path]
    //   .gwa    → DB write path, [legacy intervals path]
    //   .pheno  → DB analysis path, [legacy intervals path]
    //   .grm    → [legacy intervals path only]
    //   .plink  → [legacy intervals path only]
    // Same ConcurrentModificationException prevention as GCTA_SIMULATE_PHENOTYPES forks.
    ch_perform_gwa.params
        .tap { ch_gwa_params_for_db }
        .tap { ch_gwa_params_for_analysis }
        .set { ch_gwa_params_for_legacy }

    ch_perform_gwa.gwa
        .tap { ch_gwa_gwa_for_db }
        .set { ch_gwa_gwa_for_legacy }

    ch_perform_gwa.pheno
        .tap { ch_gwa_pheno_for_analysis }
        .set { ch_gwa_pheno_for_legacy }

    // ── DATABASE WRITE PATH (parallel fork) ──────────────────────────────
    // Writes raw GWA results to a Parquet database alongside the existing
    // QTL analysis chain. The database is a bonus artifact — if DB writes
    // fail, the pipeline still produces simulation_assessment_results.tsv.

    // ── MARKER SET CREATION ──────────────────────────────────────────
    // Wire from upstream PLINK + EIGEN channels (not from GCTA_PERFORM_GWA).
    // This eliminates the need for keyed_plink output or groupTuple dedup.
    //
    // PLINK_RECODE_MS_VCF.out.plink emits:
    //   tuple val(meta.id), val(maf), path(bed), path(bim), path(fam),
    //         path(map), path(nosex), path(ped), path(log)
    // We extract: (group, maf, bim)
    ch_bim_for_marker = PLINK_RECODE_MS_VCF.out.plink
        .map { group, maf, _bed, bim, _fam, _map_f, _nosex, _ped, _log_f ->
            tuple(group, maf, bim)
        }

    // LOCAL_COMPILE_EIGENS.out.tests emits:
    //   tuple val(group), val(maf), path(n_indep_tests)
    // ch_marker_set_params emits:
    //   [group_id, ms_maf, species, vcf_release_id, ms_ld, strains, strainfile_path]
    // Join by (group, maf) — 1:1 since all three channels emit once per key
    // Result: tuple(group, maf, bim, n_indep_tests, species, vcf_release_id, ms_ld, strains, strainfile_path)
    ch_marker_set_inputs = ch_bim_for_marker
        .join(LOCAL_COMPILE_EIGENS.out.tests, by: [0, 1])
        .join(ch_marker_set_params_for_ms, by: [0, 1])

    DB_MIGRATION_WRITE_MARKER_SET(ch_marker_set_inputs, db_output_dir)
    ch_versions = ch_versions.mix(DB_MIGRATION_WRITE_MARKER_SET.out.versions)

    // WRITE_GENOTYPE_MATRIX — join(by:[0,1]) is 1:1 per group
    // .map{} trims strains/strainfile (join artifacts not needed by write_genotype_matrix.R).
    // tuple() constructs a fresh object per emission — no shared ArrayList references (safe
    // from the ConcurrentModificationException class documented in issues #148/#154).
    ch_gm_inputs = BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix
        .join(ch_marker_set_params_for_gm, by: [0, 1])
        .map { group, maf, gm, species, vcf_release_id, ms_ld, _strains, _strainfile ->
            tuple(group, maf, gm, species, vcf_release_id, ms_ld)
        }

    DB_MIGRATION_WRITE_GENOTYPE_MATRIX(ch_gm_inputs, db_output_dir)
    ch_versions = ch_versions.mix(DB_MIGRATION_WRITE_GENOTYPE_MATRIX.out.versions)

    // ── GWA DATABASE WRITES ──────────────────────────────────────────
    // Barrier: wait for ALL marker sets to complete before writing mappings.
    //
    // .collect() on an empty channel emits [] (does NOT hang), so we
    // validate that at least one marker set was written. Without this
    // guard, a complete WRITE_MARKER_SET failure would silently allow
    // WRITE_GWA_TO_DB to proceed with no marker sets in the database.
    // ch_marker_barrier collects 2N emissions per run: N from WRITE_MARKER_SET + N from
    // WRITE_GENOTYPE_MATRIX (where N = number of population/MAF pairs in the strainfile).
    // Both signals are required so assess_sims.R can call read_genotype_matrix() safely.
    //
    // Liveness check only (size > 0), not completeness check (size == 2N).
    // If WRITE_GENOTYPE_MATRIX fails silently, ASSESS_SIMS will fail on read
    // with a descriptive tryCatch message (see assess_sims.R).
    ch_marker_barrier = DB_MIGRATION_WRITE_MARKER_SET.out.done
        .mix(DB_MIGRATION_WRITE_GENOTYPE_MATRIX.out.done)
        .collect()
        .map { items ->
            if (items.size() == 0) {
                error("No marker sets were written — cannot proceed with DB population")
            }
            true
        }

    // Gate WRITE_TRAIT_DATA behind ch_marker_barrier.
    // write_trait_data.R calls read_marker_set_metadata() — requires the marker
    // set metadata to be fully written before any trait write starts.
    // Pattern mirrors ch_db_params below: combine strips the barrier sentinel.
    ch_gated_trait_params = ch_trait.write_params
        .combine(ch_marker_barrier)
        .map { group, maf, nqtl, effect, rep, h2, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2)
        }

    DB_MIGRATION_WRITE_TRAIT_DATA(
        ch_gated_trait_params,
        ch_trait.write_pheno,
        ch_trait.write_par,
        db_output_dir,
        ch_trait.causal_geno,
        ch_trait.cv_maf,
        Channel.value(cv_ld)
    )
    ch_versions = ch_versions.mix(DB_MIGRATION_WRITE_TRAIT_DATA.out.versions)

    // Gate GWA params behind the barrier using .combine()
    //
    // .combine() with a 1-element channel preserves emission order and
    // cardinality: N params × 1 barrier = N elements in original order.
    // This delays params until the barrier resolves, while gwa (ungated)
    // waits for its corresponding params element via Nextflow's implicit
    // emission-order synchronization.
    //
    // GCTA_PERFORM_GWA.out.params emits (10 vals):
    //   tuple val(group), val(maf), val(nqtl), val(effect), val(rep),
    //         val(h2), val(mode), val(suffix), val(type), val(species)
    //
    // The DB write path intercepts GCTA_PERFORM_GWA.out BEFORE the
    // threshold expansion (Channel.of("BF", "EIGEN")). WRITE_GWA_TO_DB
    // runs once per (mode, type) combination, NOT once per threshold.
    ch_db_params = ch_gwa_params_for_db
        .combine(ch_marker_barrier)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species)
        }

    // GCTA_PERFORM_GWA.out.gwa emits a single path per invocation.
    // Nextflow pairs it with ch_gwa_db_inputs by emission index (lock-step).
    // No barrier gating needed on gwa — the params gate is sufficient.
    //
    // write_gwa_to_db.R reads species/vcf_release_id/ms_ld from marker_set_metadata.parquet
    // at runtime (via read_marker_set_metadata()). Species is still passed through the
    // channel tuple for the process tag and NF_TRAP_CELL_KEY — it's not used for ID
    // computation. The ch_marker_barrier gate on ch_db_params (above) guarantees
    // DB_MIGRATION_WRITE_MARKER_SET has completed and marker_set_metadata.parquet
    // exists before any WRITE_GWA_TO_DB task starts. If this barrier is ever weakened
    // or removed, write_gwa_to_db.R will fail with "Marker set metadata not found".
    //
    // ch_db_params: (group, maf, nqtl, effect, rep, h2, mode, suffix, type, species) — 10 elements
    // After combine(by:0) with cv_maf_keyed: adds cv_maf_eff → 11 total.
    // suffix is destructured explicitly (positional correctness) then discarded from output.
    // Species rides through from GCTA_PERFORM_GWA.out.params — no separate combine.
    ch_gwa_db_inputs = ch_db_params
        .combine(ch_cv_maf_keyed_for_gwa_write, by: 0)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species, cv_maf_eff ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, type, cv_maf_eff, cv_ld, species)
        }
    // Result: tuple(group, maf, nqtl, effect, rep, h2, mode, type, cv_maf_effective, cv_ld, species)

    DB_MIGRATION_WRITE_GWA_TO_DB(
        ch_gwa_db_inputs,
        ch_gwa_gwa_for_db,
        db_output_dir
    )
    ch_versions = ch_versions.mix(DB_MIGRATION_WRITE_GWA_TO_DB.out.versions)

    // ── METADATA AGGREGATION ─────────────────────────────────────────
    // Trigger: fire-once gate. List may be shorter than the param grid if any
    // reps were dropped under errorStrategy = 'ignore'. No size: constraint —
    // AGGREGATE_METADATA scans the filesystem directly; list contents unused.
    // Does NOT block the existing R_GET_GCTA_INTERVALS → R_ASSESS_SIMS chain.
    DB_MIGRATION_AGGREGATE_METADATA(
        DB_MIGRATION_WRITE_GWA_TO_DB.out.done.collect(),
        db_output_dir
    )
    ch_versions = ch_versions.mix(DB_MIGRATION_AGGREGATE_METADATA.out.versions)

    // ── DB-PATH QTL ANALYSIS (default) ─────────────────────────────────
    // Queries the populated database, detects QTL intervals with flexible
    // thresholds, and assesses against simulated truth. Produces
    // db_simulation_assessment_results.tsv as the primary output.

    // Barrier: wait for DB to be fully populated before querying
    ch_db_analysis_barrier = DB_MIGRATION_AGGREGATE_METADATA.out.summary

    // Expand GCTA_PERFORM_GWA params × pheno by threshold (BF × EIGEN)
    // G4: Also extend with cv_maf_effective and cv_ld so analyze_qtl.R and assess_sims.R
    // can reconstruct the v=2 trait_id (which encodes cv pool params).
    //
    // Merge params + pheno into a single channel BEFORE threshold expansion
    // so they travel together through all transforms. multiMap at the end
    // splits into synchronized sub-channels — paired by construction, not
    // emission-index coincidence. This also prevents the shared-ArrayList
    // ConcurrentModificationException from the old wrap-unwrap pattern.
    // ch_gwa_params_for_analysis carries species at position 9 (from
    // GCTA_PERFORM_GWA.out.params). No separate species combine needed —
    // species rides through forward.
    ch_db_sthresh = Channel.of("BF", "EIGEN")
    ch_gwa_params_for_analysis
        .merge(ch_gwa_pheno_for_analysis)
        // → [group, maf, nqtl, effect, rep, h2, mode, suffix, type, species, pheno, par]
        .combine(ch_db_sthresh)
        .combine(ch_db_analysis_barrier)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
               pheno, par, threshold, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                  pheno, par, threshold)
        }
        .combine(ch_cv_maf_keyed_for_analysis, by: 0)
        // → [group, maf, ..., species, pheno, par, threshold, cv_maf_eff]
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                    pheno, par, threshold, cv_maf_eff ->
            params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                          threshold, cv_maf_eff, cv_ld, species)
            pheno:  tuple(pheno, par)
        }
        .set { ch_db_analysis }

    // Step 1: Analyze QTL (pheno is pass-through for Step 2)
    DB_MIGRATION_ANALYZE_QTL(
        ch_db_analysis.params,
        ch_db_analysis.pheno,
        db_output_dir,
        params.ci_size,
        params.group_qtl,
        params.alpha
    )
    ch_versions = ch_versions.mix(DB_MIGRATION_ANALYZE_QTL.out.versions)

    // Step 2: Assess Sims (consumes ANALYZE_QTL outputs — correctly paired via pass-through)
    DB_MIGRATION_ASSESS_SIMS(
        DB_MIGRATION_ANALYZE_QTL.out.params,
        DB_MIGRATION_ANALYZE_QTL.out.regions,
        DB_MIGRATION_ANALYZE_QTL.out.pheno,
        db_output_dir,
        params.ci_size,
        params.group_qtl,
        params.alpha
    )
    ch_versions = ch_versions.mix(DB_MIGRATION_ASSESS_SIMS.out.versions)

    ch_db_assessment_pub = DB_MIGRATION_ASSESS_SIMS.out.assessment.collectFile(
        name: "db_simulation_assessment_results.tsv", sort: false
    )

    // ── REPLICATION VALIDATION ─────────────────────────────────────────
    // Runs after all DB writes have settled. Exits non-zero and writes
    // REPLAY_REQUIRED if any parameter cell is short of params.reps
    // replications or any data.parquet file is unreadable.
    def ch_validate_trigger = DB_MIGRATION_WRITE_GWA_TO_DB.out.done
        .mix(DB_MIGRATION_WRITE_TRAIT_DATA.out.done)
        .mix(DB_MIGRATION_WRITE_GENOTYPE_MATRIX.out.done)
        .mix(DB_MIGRATION_AGGREGATE_METADATA.out.summary)
        .collect()

    VALIDATE_REPLICATION_COMPLETE(
        ch_validate_trigger,
        db_output_dir,
        params.reps
    )

    // ── LEGACY ASSESSMENT PATH (optional, --legacy_assess) ────────────
    // Runs the original R_GET_GCTA_INTERVALS + R_ASSESS_SIMS chain for
    // cross-validation against the DB path. Produces
    // simulation_assessment_results.tsv alongside the DB output.
    if (params.legacy_assess) {
        // Merge all 5 upstream channels so they travel together through
        // the threshold expansion. multiMap splits into synchronized
        // sub-channels — same pattern as the DB analysis path above.
        // All originate from GCTA_PERFORM_GWA outputs (lock-step emission order).
        // After 4-alias fanout, .grm/.plink come from the multiMap-derived
        // ch_perform_gwa sub-channels — these stay aligned with the other
        // forks (.params, .pheno, .gwa) because they all emit in lockstep
        // from the same multiMap source tuple.
        // ch_gwa_params_for_legacy carries species at position 9 (from
        // GCTA_PERFORM_GWA.out.params). No separate species combine needed.
        ch_intervals_sthresh = Channel.of("BF", "EIGEN")
        ch_gwa_params_for_legacy                          // 10 vals (incl. species)
            .merge(ch_perform_gwa.grm)                    // + 3 paths (grm)
            .merge(ch_perform_gwa.plink)                  // + 9 paths (plink)
            .merge(ch_gwa_pheno_for_legacy)               // + 2 paths (pheno)
            .merge(ch_gwa_gwa_for_legacy)                 // + 1 path  (gwa)
            // → 25-element tuple
            .combine(ch_intervals_sthresh)
            // → 26 elements (+ threshold)
            .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type, species,
                        grm_bin, grm_n, grm_id,
                        bed, bim, fam, map_f, nosex, ped, log_f, gm, n_indep,
                        pheno, par,
                        gwa,
                        threshold ->
                params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, threshold, species)
                grm:    tuple(grm_bin, grm_n, grm_id)
                plink:  tuple(bed, bim, fam, map_f, nosex, ped, log_f, gm, n_indep)
                pheno:  tuple(pheno, par)
                gwa:    gwa
            }
            .set { ch_intervals }

        R_GET_GCTA_INTERVALS(
            ch_intervals.params,
            ch_intervals.grm,
            ch_intervals.plink,
            ch_intervals.pheno,
            ch_intervals.gwa,
            Channel.fromPath("${workflow.projectDir}/bin/Get_GCTA_Intervals.R").first(),
            params.group_qtl,
            params.ci_size,
            params.alpha
        )
        ch_versions = ch_versions.mix(R_GET_GCTA_INTERVALS.out.versions)

        R_ASSESS_SIMS(
            R_GET_GCTA_INTERVALS.out.params,
            R_GET_GCTA_INTERVALS.out.grm,
            R_GET_GCTA_INTERVALS.out.plink,
            R_GET_GCTA_INTERVALS.out.pheno,
            R_GET_GCTA_INTERVALS.out.interval,
            Channel.fromPath("${workflow.projectDir}/bin/Assess_Sims.R").first(),
            params.alpha,
            params.ci_size,
            params.group_qtl
        )
        ch_versions = ch_versions.mix(R_ASSESS_SIMS.out.versions)

        ch_mapping_pub = R_ASSESS_SIMS.out.assessment.collectFile(
            name: "simulation_assessment_results.tsv", sort: false
        )
    } else {
        ch_mapping_pub = Channel.empty()
    }

    // // Split results by algorithm and compile into summary file
    // R_ASSESS_SIMS.out.assessment.branch{ v ->
    //     //inbred: v.contains("inbred_nopca")
    //     inbred_pca: v.contains("inbred_pca")
    //     //loco: v.contains("loco_nopca")
    //     //loco_pca: v.contains("loco_pca")
    //     }.set{ result }

    // result.inbred_pca.view()

    publish:
        ch_mapping_pub >> "."
        ch_db_assessment_pub >> "."
        ch_input_files >> "inputs"
}


// Current bug that publish doesn't work without an output closure
output {
    "." {
        mode "copy"
    }
    "inputs" {
        mode "copy"
    }
}


/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *  GENERATE REPORT  * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/

workflow.onComplete {

    def currentSession = workflow.sessionId.toString()
    def failuresDir    = file("${workflow.outputDir}/.failures")
    def failures       = []

    if (failuresDir.exists()) {
        failuresDir.listFiles().each { f ->
            if (f.name.endsWith('.json')) {
                try {
                    def record = new groovy.json.JsonSlurper().parse(f)
                    if (record.session == currentSession) {
                        failures << record
                    }
                } catch (Exception e) {
                    log.warn "Could not parse failure record ${f.name}: ${e.message} — skipping"
                }
            }
        }
    }

    // NOTE: workflow.trace is NOT a declared property of WorkflowMetadata in
    // NF 24.10.x. Accessing it calls WorkflowMetadata.get("trace") →
    // InvokerHelper.getProperty(this,"trace") → no getTrace() getter found →
    // falls back to GroovyObject.getProperty → InvokerHelper again → infinite
    // recursion → StackOverflowError (an Error, not Exception, so it escapes
    // invokeOnComplete's catch block and crashes the session).
    //
    // The SIGKILL safety net (synthesising replay.tsv entries for walltime-killed
    // tasks) must be implemented via a TraceObserver, not workflow.trace. Until
    // that is in place, SIGKILLed tasks will not appear in replay.tsv.

    if (failures) {
        // Derive replay.tsv schema from the union of JSON keys across all failure
        // records. New fanout variables (e.g. species, future rep-disambiguators)
        // automatically appear as new columns without editing main.nf — adding a
        // field to NF_TRAP_PAYLOAD in any process flows through to replay.tsv.
        //
        // Column ordering: stable across runs. Fields present in any record but
        // missing from a particular record render as empty cells.
        def columnOrder = [
            'session', 'task_hash', 'attempt', 'max_retries',
            'species', 'group', 'maf', 'nqtl', 'effect', 'h2', 'rep', 'mode', 'type',
            'exit'
        ]
        def discoveredKeys = failures.collectMany { it.keySet() as List }.unique()
        def orderedKnown   = columnOrder.findAll { it in discoveredKeys }
        def extras         = (discoveredKeys - columnOrder).sort()
        def headerCols     = orderedKnown + extras

        def lines = new StringBuilder()
        lines << headerCols.join('\t') << '\n'
        failures.each { r ->
            lines << headerCols.collect { col -> r.containsKey(col) ? r[col] : '' }.join('\t') << '\n'
        }
        def tsvTmp = file("${workflow.outputDir}/replay.tsv.tmp")
        def tsv    = file("${workflow.outputDir}/replay.tsv")
        tsvTmp.text = lines.toString()
        tsvTmp.renameTo(tsv)

        def byExit   = failures.groupBy { it.exit }
        def exitSummary = byExit.collect { exit, list -> "${list.size()} × exit ${exit}" }.join(', ')
        log.warn "Replay manifest written: ${failures.size()} failed slots (${exitSummary}) → ${tsv}"
    } else {
        log.info "No failures recorded — replay.tsv not written."
    }

    // Pre-capture git properties — accessing workflow.repository/revision/commitId
    // inside a GString via $workflow.xxx triggers infinite recursion in
    // WorkflowMetadata.get() in NF 24.10.x (Groovy dynamic dispatch loop).
    // Fall back to local `git` commands for local runs (workflow fields are null
    // when the pipeline is launched from a directory rather than a GitHub URL).
    def gitRepo     = workflow.repository ?: ("git remote get-url origin".execute().text.trim() ?: 'N/A')
    def gitRevision = workflow.revision   ?: ("git rev-parse --abbrev-ref HEAD".execute().text.trim() ?: 'N/A')
    def gitCommit   = workflow.commitId   ?: ("git rev-parse HEAD".execute().text.trim() ?: 'N/A')

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: ${gitRepo} - ${gitRevision} [${gitCommit}]

    { Parameters }
    ---------------------------
    Strainfile                              = ${params.strainfile}
    Causal variant MAF                      = ${cv_maf ?: '(per-group ms_maf)'}
    Causal variant LD threshold             = ${cv_ld}
    Number of simulated QTLs                = ${nqtl_file}
    Phenotype Heritability File             = ${h2_file}
    Number of simulation replicates         = ${params.reps}
    Effect Size Range File                  = ${effect_file}
    Marker Genomic Range File               = ${params.qtlloc}
    Significance Thresholds                 = BF, EIGEN
    Threshold for grouping QTL              = ${params.group_qtl}
    Number of SNVs to define CI             = ${params.ci_size}
    Relatedness Matrix Cutoff               = ${params.sparse_cut}
    Mitochondrial chromosome name           = ${mito_name}
    Simulate QTLs in specific regions       = ${simulate_qtlloc}
    Result Directory                        = ${workflow.outputDir}
    """

    println summary

}
