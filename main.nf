// Needed to publish results (NF 24.10.x preview feature)
nextflow.preview.output = true

// import the subworkflows
include { LOCAL_GET_CONTIG_INFO           } from './modules/local/get_contig_info/main'
include { LOCAL_COMPILE_EIGENS            } from './modules/local/compile_eigens/main'
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
include { GCTA_MAKE_GRM                   } from './modules/gcta/make_grm/main'
include { GCTA_PERFORM_GWA                } from './modules/gcta/perform_gwa/main'

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
        switch (species) {
            case 'c_elegans':    return "https://caendr.org/download/WI.${vcf}.hard-filter.isotype.vcf.gz"
            case 'c_briggsae':   return "https://caendr.org/download/CB.${vcf}.hard-filter.isotype.vcf.gz"
            case 'c_tropicalis': return "https://caendr.org/download/CT.${vcf}.hard-filter.isotype.vcf.gz"
            default:
                throw new IllegalArgumentException(
                    "Unrecognized species '${species}' for CaeNDR URL construction. " +
                    "Use 'c_elegans', 'c_briggsae', or 'c_tropicalis'."
                )
        }
    }
    if (!vcf.startsWith('/') && !vcf.startsWith('http')) {
        return "${projectDir}/${vcf}"  // relative path → absolute (test profile support)
    }
    return vcf  // absolute path passthrough
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
        log.info "Mitochondrial chromosome name           = ${mito_name}"
        log.info "Simulate QTLs in specific regions       = ${simulate_qtlloc}"
        log.info "Output directory                        = ${workflow.outputDir}"
        log.info ""
    }



    // Parse 6-column strainfile and fan out per-row parameters
    ch_strain_sets = Channel.fromPath(strainfile)
        .splitCsv(sep: "\t", header: true)
        .map { row ->
            // M8 fix: normalize whitespace to handle Windows CRLF line endings
            def cleanRow = row.collectEntries { k, v -> [k.trim(), v?.trim()] }
            if (!(cleanRow.species in ['c_elegans', 'c_briggsae', 'c_tropicalis'])) {
                error "Unknown species '${cleanRow.species}' in strainfile row for group '${cleanRow.group}'"
            }
            if (cleanRow.ms_maf.toFloat() <= 0 || cleanRow.ms_maf.toFloat() > 0.5) {
                error "ms_maf '${cleanRow.ms_maf}' out of range (0, 0.5] in group '${cleanRow.group}' — values above 0.5 are major allele frequency, not MAF"
            }
            if (cleanRow.ms_ld.toFloat() <= 0 || cleanRow.ms_ld.toFloat() >= 1.0) {
                error "ms_ld '${cleanRow.ms_ld}' out of range (0, 1) in group '${cleanRow.group}'"
            }
            [
                [id: cleanRow.group],
                cleanRow.species,
                cleanRow.vcf,
                cleanRow.ms_maf.toFloat(),
                cleanRow.ms_ld.toFloat(),
                cleanRow.strains
            ]
        }

    ch_strain_sets
        .multiMap { meta, species, vcf, ms_maf, ms_ld, strains ->
            def cv_maf_eff = cv_maf != null ? cv_maf : ms_maf
            if (cv_maf != null && cv_maf > ms_maf) {
                log.warn "cv_maf (${cv_maf}) > ms_maf (${ms_maf}) for group ${meta.id}: " +
                         "the CV pool is a strict subset of the marker SNP set. " +
                         "This is likely unintentional."
            }
            marker_set_params: [meta.id, ms_maf, species, extractVcfReleaseId(vcf), ms_ld, strains, strainfile]
            vcf_per_group:     [meta, species, vcf, strains]
            ms_maf_vals:       [meta, ms_maf]
            ms_ld_vals:        [meta, ms_ld]
            cv_maf_vals:       [meta, cv_maf_eff]
        }
        .set { ch_sf }

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
        .map { meta, species, vcf, strains ->
            def vcf_path = resolveVcf(vcf, species, workflow.projectDir.toString())
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
    PYTHON_SIMULATE_EFFECTS_GLOBAL(
        ch_plink_genomat_eigen,
        Channel.fromPath("${workflow.projectDir}/bin/create_causal_vars.py").first(),
        Channel.of(1..params.reps).toSortedList(),
        Channel.fromPath(nqtl_file)
            .splitCsv()
            .map{ it: it[0] }
            .toSortedList(),
        Channel.fromPath(effect_file)
            .splitCsv()
            .map{ it: it[0] }
            .toSortedList() 
        )
    ch_versions = ch_versions.mix(PYTHON_SIMULATE_EFFECTS_GLOBAL.out.versions)
    ch_sim_phenos    = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.causal
    ch_sim_plink     = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.plink
    ch_sim_cv_plink  = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.cv_plink
    // }

    // G2: Fan causal_genotypes across h2 to match GCTA_SIMULATE_PHENOTYPES cardinality.
    // causal_genotypes emits once per (group, maf, nqtl, effect, rep) — before h2 expansion.
    // Each geno file must be paired with every h2 value so the merge chain join (G3) can
    // align per (group, maf, nqtl, effect, rep, h2).
    ch_causal_geno_fanned = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.causal_genotypes
        .combine(Channel.fromPath(h2_file).splitCsv().map { it[0] })
        .map { group, maf, nqtl, effect, rep, geno_file, h2 ->
            tuple(group, maf, nqtl, effect, rep, h2, geno_file)
        }

    // Simulate phenotypes using the CV binary (--bfile CV_TO_SIMS) so GCTA can locate
    // non-marker causal variant genotypes. The MS binary (ch_sim_plink) passes through
    // to downstream GWA mapping unchanged.
    GCTA_SIMULATE_PHENOTYPES(
        ch_sim_phenos,
        ch_sim_cv_plink,   // CV binary — used as --bfile for GCTA simulation
        ch_sim_plink,      // MS binary — passed through to downstream GWA
        Channel.fromPath("${h2_file}")
            .splitCsv()
            .map{ it: it[0] }
            .toSortedList()
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

    // Resolve db_output to absolute path so SLURM tasks write to the correct
    // shared filesystem location, not relative to their work directory.
    // Default: {outputDir}/db (computed from workflow.outputDir when params.db_output is null)
    def db_output_dir = params.db_output
        ? file(params.db_output).toAbsolutePath().toString()
        : file("${workflow.outputDir}/db").toAbsolutePath().toString()
    log.info "Database output directory: ${db_output_dir}"

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
    // After .combine(ch_cv_maf_keyed_for_trait, by: 0): 20-element tuple (cv_maf_eff at [19])
    // After .join(ch_causal_geno_fanned, by: [0..5]):   21-element tuple (geno_file at [20])
    //
    // WARNING: .merge() aligns channels by emission order, NOT by key matching.
    // Do NOT insert .filter(), .map(), .branch(), or any reordering operator
    // between GCTA_SIMULATE_PHENOTYPES.out.* and this .merge() chain — doing
    // so will silently misalign phenotype/causal files with metadata.
    // Runtime validation in write_trait_data.R catches misalignment.
    // .combine(by:0) and .join(by:[0..5]) are key-based (order-preserving) and
    // therefore safe per the CLAUDE.md merge ordering constraint.
    //
    // Merged tuple structure (21 elements):
    //   [0-5]   params:   group, maf, nqtl, effect, rep, h2
    //   [6-7]   pheno:    phen_path, par_path
    //   [8-18]  plink:    group, maf, bed, bim, fam, map, nosex, ped, log, gm, n_indep_tests
    //   [19]    cv_maf:   cv_maf_effective (from combine)
    //   [20]    geno:     causal_genotypes TSV path (from join)
    //
    // IMPORTANT: If GCTA_SIMULATE_PHENOTYPES output channels change, update
    // these indices.
    ch_gcta_params_for_trait
        .merge(ch_gcta_pheno_for_trait)
        .merge(ch_gcta_plink_for_trait)
        .combine(ch_cv_maf_keyed_for_trait, by: 0)
        .join(ch_causal_geno_fanned, by: [0, 1, 2, 3, 4, 5])
        .multiMap { it ->
            assert it.size() == 21 : "Expected 21-element tuple, got ${it.size()}. First elements: ${it.take(6)}. Check output channels in modules/gcta/simulate_phenotypes/main.nf"
            // trait_id is computed in R by write_trait_data.R via generate_trait_id()
            // — no cross-language hash computation (addresses review B2)
            write_params:  tuple(it[0], it[1], it[2], it[3], it[4], it[5])
            write_pheno:   it[6]   // pre-upscaled .phen file
            write_par:     it[7]   // causal variant .par file
            cv_maf:        it[19]  // cv_maf_effective (from combine with ch_cv_maf_keyed_for_trait)
            causal_geno:   it[20]  // causal genotype TSV (from join with ch_causal_geno_fanned)
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

    // Create genetic relatedness matrix
    // Merge PLINK_UPDATE_BY_H2 outputs before the mode fan-out (inbred/loco).
    // The previous wrap-combine-unwrap pattern (.map{[it]}.combine(ch_mode).map{it[0]})
    // returned the same ArrayList reference for both mode emissions. Under SLURM job
    // array batching + Singularity, two submission threads iterated the shared list
    // simultaneously -> ConcurrentModificationException (same mechanism as issue #148).
    // .merge() aligns by emission order — do not insert reordering operators between
    // PLINK_UPDATE_BY_H2.out.* and this merge chain.
    ch_mode = Channel.of(
        ["inbred", "fastGWA"],
        ["loco", "mlma"]
        )
    PLINK_UPDATE_BY_H2.out.params
        .merge(PLINK_UPDATE_BY_H2.out.plink)
        .merge(PLINK_UPDATE_BY_H2.out.pheno)
        .combine(ch_mode)
        .multiMap { group, maf, nqtl, effect, rep, h2,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par, mode, suffix ->
            params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix)
            plink:  tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno:  tuple(pheno, par)
        }
        .set { ch_grm }

    GCTA_MAKE_GRM(
        ch_grm.params,
        ch_grm.plink,
        ch_grm.pheno
        )
    ch_versions = ch_versions.mix(GCTA_MAKE_GRM.out.versions)

    // Simulate GWA using output from GCTA_MAKE_GRM
    // Same fix: merge all GCTA_MAKE_GRM outputs before the type fan-out (pca/nopca).
    // .merge() aligns by emission order — do not insert reordering operators between
    // GCTA_MAKE_GRM.out.* and this merge chain.
    ch_type = Channel.of(
        "pca",
        "nopca"
        )
    GCTA_MAKE_GRM.out.params
        .merge(GCTA_MAKE_GRM.out.grm)
        .merge(GCTA_MAKE_GRM.out.plink)
        .merge(GCTA_MAKE_GRM.out.pheno)
        .combine(ch_type)
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix,
                    grm_bin, grm_n, grm_id,
                    bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                    pheno, par, type ->
            params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type)
            grm:    tuple(grm_bin, grm_n, grm_id)
            plink:  tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
            pheno:  tuple(pheno, par)
        }
        .set { ch_gwa }

    GCTA_PERFORM_GWA(
        ch_gwa.params,
        ch_gwa.grm,
        ch_gwa.plink,
        ch_gwa.pheno
        )
    ch_versions = ch_versions.mix(GCTA_PERFORM_GWA.out.versions)

    // Fork GCTA_PERFORM_GWA outputs for multiple consumers:
    //   .out.params → DB write path, DB analysis path, [legacy intervals path]
    //   .out.gwa    → DB write path, [legacy intervals path]
    //   .out.pheno  → DB analysis path, [legacy intervals path]
    // Same ConcurrentModificationException prevention as GCTA_SIMULATE_PHENOTYPES forks.
    GCTA_PERFORM_GWA.out.params
        .tap { ch_gwa_params_for_db }
        .tap { ch_gwa_params_for_analysis }
        .set { ch_gwa_params_for_legacy }

    GCTA_PERFORM_GWA.out.gwa
        .tap { ch_gwa_gwa_for_db }
        .set { ch_gwa_gwa_for_legacy }

    GCTA_PERFORM_GWA.out.pheno
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
    ch_gm_inputs = BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix
        .join(ch_marker_set_params_for_gm, by: [0, 1])
    // Result: tuple(group, maf, genotype_matrix, species, vcf_release_id, ms_ld, strains, strainfile)

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
    // GCTA_PERFORM_GWA.out.params emits:
    //   tuple val(group), val(maf), val(nqtl), val(effect), val(rep),
    //         val(h2), val(mode), val(suffix), val(type)
    //
    // The DB write path intercepts GCTA_PERFORM_GWA.out BEFORE the
    // threshold expansion (Channel.of("BF", "EIGEN")). WRITE_GWA_TO_DB
    // runs once per (mode, type) combination, NOT once per threshold.
    ch_db_params = ch_gwa_params_for_db
        .combine(ch_marker_barrier)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type)
        }

    // GCTA_PERFORM_GWA.out.gwa emits a single path per invocation.
    // Nextflow pairs it with ch_gwa_db_inputs by emission index (lock-step).
    // No barrier gating needed on gwa — the params gate is sufficient.
    //
    // write_gwa_to_db.R reads species/vcf_release_id/ms_ld from marker_set_metadata.parquet
    // at runtime (via read_marker_set_metadata()). The ch_marker_barrier gate on ch_db_params
    // (above) guarantees DB_MIGRATION_WRITE_MARKER_SET has completed and
    // marker_set_metadata.parquet exists before any WRITE_GWA_TO_DB task starts.
    // If this barrier is ever weakened or removed, write_gwa_to_db.R will fail
    // with "Marker set metadata not found" errors.
    //
    // ch_db_params: (group, maf, nqtl, effect, rep, h2, mode, suffix, type) — 9 elements
    // After combine(by:0) with cv_maf_keyed: adds cv_maf_eff → 10 total.
    // suffix is destructured explicitly (positional correctness) then discarded from output.
    ch_gwa_db_inputs = ch_db_params
        .combine(ch_cv_maf_keyed_for_gwa_write, by: 0)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, cv_maf_eff ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, type, cv_maf_eff, cv_ld)
        }
    // Result: tuple(group, maf, nqtl, effect, rep, h2, mode, type, cv_maf_effective, cv_ld)

    DB_MIGRATION_WRITE_GWA_TO_DB(
        ch_gwa_db_inputs,
        ch_gwa_gwa_for_db,
        db_output_dir
    )
    ch_versions = ch_versions.mix(DB_MIGRATION_WRITE_GWA_TO_DB.out.versions)

    // ── METADATA AGGREGATION ─────────────────────────────────────────
    // Runs after ALL WRITE_GWA_TO_DB processes complete.
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
    ch_db_sthresh = Channel.of("BF", "EIGEN")
    ch_gwa_params_for_analysis
        .merge(ch_gwa_pheno_for_analysis)
        // → [group, maf, nqtl, effect, rep, h2, mode, suffix, type, pheno, par]
        .combine(ch_db_sthresh)
        .combine(ch_db_analysis_barrier)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
               pheno, par, threshold, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                  pheno, par, threshold)
        }
        .combine(ch_cv_maf_keyed_for_analysis, by: 0)
        // → [group, maf, ..., pheno, par, threshold, cv_maf_eff]
        .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                    pheno, par, threshold, cv_maf_eff ->
            params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                          threshold, cv_maf_eff, cv_ld)
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

    // ── LEGACY ASSESSMENT PATH (optional, --legacy_assess) ────────────
    // Runs the original R_GET_GCTA_INTERVALS + R_ASSESS_SIMS chain for
    // cross-validation against the DB path. Produces
    // simulation_assessment_results.tsv alongside the DB output.
    if (params.legacy_assess) {
        // Merge all 5 upstream channels so they travel together through
        // the threshold expansion. multiMap splits into synchronized
        // sub-channels — same pattern as the DB analysis path above.
        // All originate from GCTA_PERFORM_GWA outputs (lock-step emission order).
        ch_intervals_sthresh = Channel.of("BF", "EIGEN")
        ch_gwa_params_for_legacy                          // 9 vals
            .merge(GCTA_PERFORM_GWA.out.grm)              // + 3 paths (grm)
            .merge(GCTA_PERFORM_GWA.out.plink)            // + 9 paths (plink)
            .merge(ch_gwa_pheno_for_legacy)               // + 2 paths (pheno)
            .merge(ch_gwa_gwa_for_legacy)                 // + 1 path  (gwa)
            // → 24-element tuple
            .combine(ch_intervals_sthresh)
            // → 25 elements (+ threshold)
            .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                        grm_bin, grm_n, grm_id,
                        bed, bim, fam, map_f, nosex, ped, log_f, gm, n_indep,
                        pheno, par,
                        gwa,
                        threshold ->
                params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, threshold)
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

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

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
    Mitochondrial chromosome name           = ${mito_name}
    Simulate QTLs in specific regions       = ${simulate_qtlloc}
    Result Directory                        = ${workflow.outputDir}
    """

    println summary

}
