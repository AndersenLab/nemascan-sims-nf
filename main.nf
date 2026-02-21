// Needed to publish results (NF 24.10.x preview feature)
nextflow.preview.output = true

// import the subworkflows
include { LOCAL_GET_CONTIG_INFO           } from './modules/local/get_contig_info/main'
include { LOCAL_COMPILE_EIGENS            } from './modules/local/compile_eigens/main'
include { BCFTOOLS_RENAME_CHROMS          } from './modules/bcftools/rename_chroms/main'
include { BCFTOOLS_EXTRACT_STRAINS        } from './modules/bcftools/extract_strains/main'
include { BCFTOOLS_CREATE_GENOTYPE_MATRIX } from './modules/bcftools/create_genotype_matrix/main'
include { PLINK_RECODE_VCF                } from './modules/plink/recode_vcf/main'
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
include { PYTHON_CHECK_VP               } from './modules/python/check_vp.nf'

// Database migration modules
include { DB_MIGRATION_WRITE_MARKER_SET       } from './modules/db_migration/write_marker_set/main'
include { DB_MIGRATION_WRITE_GWA_TO_DB        } from './modules/db_migration/write_gwa_to_db/main'
include { DB_MIGRATION_AGGREGATE_METADATA     } from './modules/db_migration/aggregate_metadata/main'
include { DB_MIGRATION_ANALYZE_QTL            } from './modules/db_migration/analyze_qtl/main'
include { DB_MIGRATION_ASSESS_SIMS            } from './modules/db_migration/assess_sims/main'
include { DB_MIGRATION_WRITE_GENOTYPE_MATRIX  } from './modules/db_migration/write_genotype_matrix/main'
include { DB_MIGRATION_WRITE_TRAIT_DATA       } from './modules/db_migration/write_trait_data/main'

workflow {
    main:
    ch_versions = Channel.empty()

    date = new Date().format( 'yyyyMMdd' )


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
    if (params.maf == null){
        maf_file = "${workflow.projectDir}/data/simulate_maf.csv"
    } else {
        maf_file = params.maf
    }
    if (params.effect == null){
        effect_file = "${workflow.projectDir}/data/simulate_effect_sizes.csv"
    } else {
        effect_file = params.effect
    }
    if (params.strainfile == null){
        strainfile = "${workflow.projectDir}/data/test_strains.txt"
    } else {
        strainfile = params.strainfile
    }
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
        log.info "--strainfile      File               A TSV file with two columns: the first is a name for the strain set and the second is a comma-separated strain list without spaces"
        log.info "--vcf             File               Generally a CaeNDR release date (i.e. 20231213). Can also provide a user-specified VCF with index in same folder"
        log.info " "
        log.info "Optional arguments (General):"
        log.info "--nqtl            File               A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: data/simulate_nqtl.csv)"
        log.info "--h2              File               A CSV file with phenotype heritability, one value per line (Default is located: data/simulate_h2.csv)"
        log.info "--reps             Integer            The number of replicates to simulate per number of QTL and heritability (Default: 2)"
        log.info "--maf             File               A CSV file where each line is a minor allele frequency threshold to test for simulations (Default: data/simulate_maf.csv)"
        log.info "--effect          File               A CSV file where each line is an effect size range (e.g. 0.2-0.3) to test for simulations (Default: data/simulate_effect_sizes.csv)"
        log.info "--qtlloc          File               A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])"
        log.info "--sthresh         String             Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
        log.info "--group_qtl       Integer            If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
        log.info "--ci_size         Integer            Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
        log.info "--sparse_cut      Decimal            Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
        log.info "--mito_name       Strain             Name of mitochondrial chromosome"
        log.info "--simulate_qtlloc Boolean            Whether to simulate QTLs in specific genomic regions (Default: false)"
        log.info "-output-dir       String             Name of folder that will contain the results (Default: Analysis_Results_{date})"
        log.info " "


        exit 1
    } else { // set log info
        log.info '''


    '''
        log.info ""
        log.info "Strain name and list file               = ${strainfile}"
        log.info "VCF                                     = ${params.vcf}"
        log.info "Number of QTLs/phenotype simulated      = ${nqtl_file}"
        log.info "Phenotype heritability file             = ${h2_file}"
        log.info "Number of replicates to simulate        = ${params.reps}"
        log.info "Minor allele freq. threshold file       = ${maf_file}"
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



    // Created needed channels
    ch_vcf = Channel.fromPath(params.vcf).map{ it: [[id: "vcf"], it, "${it}.tbi"] }.first()
    ch_strain_sets = Channel.fromPath(strainfile).splitCsv(sep: " ").map{ it: [[id: it[0]], it[1]] }
    ch_mafs = Channel.fromPath(maf_file).splitCsv().first()

    // Get contig data from VCF file
    LOCAL_GET_CONTIG_INFO( ch_vcf )
    ch_mito_num = LOCAL_GET_CONTIG_INFO.out.mapping.splitCsv(sep:"\t")
        .filter{ row -> row[0] == mito_name }
        .map{ row -> row[1] }
        .first()

    BCFTOOLS_EXTRACT_STRAINS( ch_vcf, ch_strain_sets )
    ch_versions = ch_versions.mix(BCFTOOLS_EXTRACT_STRAINS.out.versions)

    // Extract desired strain sets
    BCFTOOLS_RENAME_CHROMS( 
        BCFTOOLS_EXTRACT_STRAINS.out.vcf,
        LOCAL_GET_CONTIG_INFO.out.mapping
        )
    ch_versions = ch_versions.mix(BCFTOOLS_RENAME_CHROMS.out.versions)

    // Recode the VCF file and create plink formatted files
    PLINK_RECODE_VCF(
        BCFTOOLS_RENAME_CHROMS.out.vcf,
        ch_mito_num,
        ch_mafs
        )
    ch_versions = ch_versions.mix(PLINK_RECODE_VCF.out.versions)

    // Create plaintext genotype matrix
    BCFTOOLS_CREATE_GENOTYPE_MATRIX(
        PLINK_RECODE_VCF.out.vcf,
        PLINK_RECODE_VCF.out.markers
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

    // Compile required files for simulations by strain group and MAF
    ch_plink_genomat_eigen = PLINK_RECODE_VCF.out.plink
        .join(BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix, by: [0, 1])
        //.join(ch_eigens, by: [0, 1])
        .join(LOCAL_COMPILE_EIGENS.out.tests, by: [0, 1])
    
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
    ch_sim_phenos = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.causal
    ch_sim_plink = PYTHON_SIMULATE_EFFECTS_GLOBAL.out.plink
    // }

    // Adjust plink data by heritability
    GCTA_SIMULATE_PHENOTYPES(
        ch_sim_phenos,
        ch_sim_plink,
        Channel.fromPath("${h2_file}")
            .splitCsv()
            .map{ it: it[0] }
            .toSortedList()
        )
    ch_versions = ch_versions.mix(GCTA_SIMULATE_PHENOTYPES.out.versions)

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
    // Trait data (metadata, causal variants, phenotype) and genotype matrix are
    // stored in the DB for future analysis. These data are NOT consumed by
    // ASSESS_SIMS -- no barrier is needed.
    //
    // WARNING: .merge() aligns channels by emission order, NOT by key matching.
    // Do NOT insert .filter(), .branch(), .map(), or any reordering operator
    // between GCTA_SIMULATE_PHENOTYPES.out.* and this .merge() chain — doing
    // so will silently misalign phenotype/causal files with metadata.
    // Runtime validation in write_trait_data.R catches misalignment (see M4).
    //
    // Merged tuple structure (19 elements):
    //   [0-5]   params:  group, maf, nqtl, effect, rep, h2
    //   [6-7]   pheno:   phen_path, par_path
    //   [8-18]  plink:   group, maf, bed, bim, fam, map, nosex, ped, log, gm, n_indep_tests
    //
    // IMPORTANT: If GCTA_SIMULATE_PHENOTYPES output channels change, update
    // these indices. See "Tuple index verification" in step7-plan.qmd.

    GCTA_SIMULATE_PHENOTYPES.out.params
        .merge(GCTA_SIMULATE_PHENOTYPES.out.pheno)
        .merge(GCTA_SIMULATE_PHENOTYPES.out.plink)
        .multiMap { it ->
            assert it.size() == 19 : "Expected 19-element tuple, got ${it.size()}. First elements: ${it.take(6)}. Check output channels in modules/gcta/simulate_phenotypes/main.nf"
            // trait_id is computed in R by write_trait_data.R via generate_trait_id()
            // — no cross-language hash computation (addresses review B2)
            write_params:  tuple(it[0], it[1], it[2], it[3], it[4], it[5])
            write_pheno:   it[6]   // pre-upscaled .phen file
            write_par:     it[7]   // causal variant .par file
        }
        .set { ch_trait }

    // Write trait metadata + causal variants + phenotype to DB
    DB_MIGRATION_WRITE_TRAIT_DATA(
        ch_trait.write_params,
        ch_trait.write_pheno,
        ch_trait.write_par,
        db_output_dir
    )

    // Update plink data by heritability
    PLINK_UPDATE_BY_H2(
        GCTA_SIMULATE_PHENOTYPES.out.params,
        GCTA_SIMULATE_PHENOTYPES.out.plink,
        GCTA_SIMULATE_PHENOTYPES.out.pheno
        )
    ch_versions = ch_versions.mix(PLINK_UPDATE_BY_H2.out.versions)

    // Create genetic relatedness matrix
    ch_mode = Channel.of( 
        ["inbred", "fastGWA"],
        ["loco", "mlma"]
        )
    ch_grm_params = PLINK_UPDATE_BY_H2.out.params.combine(ch_mode)
    ch_grm_plink = PLINK_UPDATE_BY_H2.out.plink.map{ it: [it] }.combine(ch_mode).map{ it: it[0] }
    ch_grm_pheno = PLINK_UPDATE_BY_H2.out.pheno.map{ it: [it] }.combine(ch_mode).map{ it: it[0] }

    
    GCTA_MAKE_GRM(
        ch_grm_params,
        ch_grm_plink,
        ch_grm_pheno
        )
    ch_versions = ch_versions.mix(GCTA_MAKE_GRM.out.versions)

    // Prepare inputs for PYTHON_CHECK_VP to match its 3-input definition
    // Input 1: meta_tuple (group, maf, nqtl, ...)
    // Input 2: tuple (tmp_pheno_path, hsq_path, par_path)
    // Input 3: script_path

    PYTHON_CHECK_VP(
        GCTA_MAKE_GRM.out.params,                         // Arg 1: Metadata
        GCTA_MAKE_GRM.out.pheno_hsq_and_par,              // Arg 2: Tuple of (tmp_pheno, hsq, par)
        GCTA_MAKE_GRM.out.grm,
        GCTA_MAKE_GRM.out.plink,
        Channel.fromPath("${workflow.projectDir}/bin/check_vp.py").first() // Arg 3: Script path
        )
    ch_versions = ch_versions.mix(PYTHON_CHECK_VP.out.versions)

    // Simulate GWA using output from PYTHON_CHECK_VP
    ch_type = Channel.of( 
        "pca", //just PCA for now
        "nopca"
        )
    
    // Pheno file for GWA now comes from PYTHON_CHECK_VP.out.pheno
    // GCTA_PERFORM_GWA expects a tuple: (phen_path, par_path).
    // PYTHON_CHECK_VP.out.pheno is: tuple (final_pheno_path, par_path)
    // The rest of the outputs (GRM, plink, and params) are passed through
    // PYTHON_CHECK_VP to maintain correct ordering of parallel channels
    ch_gwa_params = PYTHON_CHECK_VP.out.params.combine(ch_type)
    ch_gwa_grm =    PYTHON_CHECK_VP.out.grm.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    ch_gwa_plink =  PYTHON_CHECK_VP.out.plink.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    ch_gwa_pheno =  PYTHON_CHECK_VP.out.pheno.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    
    GCTA_PERFORM_GWA(
        ch_gwa_params,
        ch_gwa_grm,
        ch_gwa_plink,
        ch_gwa_pheno,
        params.sparse_cut
        )
    ch_versions = ch_versions.mix(GCTA_PERFORM_GWA.out.versions)

    // ── DATABASE WRITE PATH (parallel fork) ──────────────────────────────
    // Writes raw GWA results to a Parquet database alongside the existing
    // QTL analysis chain. The database is a bonus artifact — if DB writes
    // fail, the pipeline still produces simulation_assessment_results.tsv.

    // ── MARKER SET CREATION ──────────────────────────────────────────
    // Wire from upstream PLINK + EIGEN channels (not from GCTA_PERFORM_GWA).
    // This eliminates the need for keyed_plink output or groupTuple dedup.
    //
    // PLINK_RECODE_VCF.out.plink emits:
    //   tuple val(meta.id), val(maf), path(bed), path(bim), path(fam),
    //         path(map), path(nosex), path(ped), path(log)
    // We extract: (group, maf, bim)
    ch_bim_for_marker = PLINK_RECODE_VCF.out.plink
        .map { group, maf, _bed, bim, _fam, _map_f, _nosex, _ped, _log_f ->
            tuple(group, maf, bim)
        }

    // LOCAL_COMPILE_EIGENS.out.tests emits:
    //   tuple val(group), val(maf), path(n_indep_tests)
    // Join by (group, maf) — 1:1 since both channels emit once per key
    ch_marker_set_inputs = ch_bim_for_marker
        .join(LOCAL_COMPILE_EIGENS.out.tests, by: [0, 1])
    // Result: tuple(group, maf, bim, n_indep_tests)

    DB_MIGRATION_WRITE_MARKER_SET(ch_marker_set_inputs, db_output_dir)

    DB_MIGRATION_WRITE_GENOTYPE_MATRIX(
        BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix,
        db_output_dir
    )

    // ── GWA DATABASE WRITES ──────────────────────────────────────────
    // Barrier: wait for ALL marker sets to complete before writing mappings.
    //
    // .collect() on an empty channel emits [] (does NOT hang), so we
    // validate that at least one marker set was written. Without this
    // guard, a complete WRITE_MARKER_SET failure would silently allow
    // WRITE_GWA_TO_DB to proceed with no marker sets in the database.
    ch_marker_barrier = DB_MIGRATION_WRITE_MARKER_SET.out.done
        .collect()
        .map { items ->
            if (items.size() == 0) {
                error("No marker sets were written — cannot proceed with DB population")
            }
            true
        }

    // Gate params behind the barrier using .combine()
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
    ch_db_params = GCTA_PERFORM_GWA.out.params
        .combine(ch_marker_barrier)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type)
        }

    // GCTA_PERFORM_GWA.out.gwa emits a single path per invocation.
    // Nextflow pairs it with ch_db_params by emission index (lock-step).
    // No barrier gating needed on gwa — the params gate is sufficient.
    DB_MIGRATION_WRITE_GWA_TO_DB(
        ch_db_params,
        GCTA_PERFORM_GWA.out.gwa,
        db_output_dir
    )

    // ── METADATA AGGREGATION ─────────────────────────────────────────
    // Runs after ALL WRITE_GWA_TO_DB processes complete.
    // Does NOT block the existing R_GET_GCTA_INTERVALS → R_ASSESS_SIMS chain.
    DB_MIGRATION_AGGREGATE_METADATA(
        DB_MIGRATION_WRITE_GWA_TO_DB.out.done.collect(),
        db_output_dir
    )

    // ── DB-PATH QTL ANALYSIS (default) ─────────────────────────────────
    // Queries the populated database, detects QTL intervals with flexible
    // thresholds, and assesses against simulated truth. Produces
    // db_simulation_assessment_results.tsv as the primary output.

    // Barrier: wait for DB to be fully populated before querying
    ch_db_analysis_barrier = DB_MIGRATION_AGGREGATE_METADATA.out.summary

    // Expand GCTA_PERFORM_GWA params by threshold (BF × EIGEN)
    ch_db_sthresh = Channel.of("BF", "EIGEN")
    ch_db_analysis_params = GCTA_PERFORM_GWA.out.params
        .combine(ch_db_sthresh)
        .combine(ch_db_analysis_barrier)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type, threshold, _barrier ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, threshold)
        }

    // Expand pheno (contains .par file) by threshold — same pattern
    ch_db_analysis_pheno = GCTA_PERFORM_GWA.out.pheno
        .map { it -> [it] }
        .combine(ch_db_sthresh)
        .combine(ch_db_analysis_barrier)
        .map { it -> it[0] }

    // Step 1: Analyze QTL (pheno is pass-through for Step 2)
    DB_MIGRATION_ANALYZE_QTL(
        ch_db_analysis_params,
        ch_db_analysis_pheno,
        db_output_dir,
        params.ci_size,
        params.group_qtl,
        params.alpha
    )

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

    ch_db_assessment_pub = DB_MIGRATION_ASSESS_SIMS.out.assessment.collectFile(
        name: "db_simulation_assessment_results.tsv", sort: false
    )

    // ── LEGACY ASSESSMENT PATH (optional, --legacy_assess) ────────────
    // Runs the original R_GET_GCTA_INTERVALS + R_ASSESS_SIMS chain for
    // cross-validation against the DB path. Produces
    // simulation_assessment_results.tsv alongside the DB output.
    if (params.legacy_assess) {
        ch_intervals_sthresh = Channel.of("BF", "EIGEN")
        ch_intervals_params = GCTA_PERFORM_GWA.out.params.combine(ch_intervals_sthresh)
        ch_intervals_grm = GCTA_PERFORM_GWA.out.grm.map{ it: [it] }.combine(ch_intervals_sthresh).map{ it: it[0] }
        ch_intervals_plink = GCTA_PERFORM_GWA.out.plink.map{ it: [it] }.combine(ch_intervals_sthresh).map{ it: it[0] }
        ch_intervals_pheno = GCTA_PERFORM_GWA.out.pheno.map{ it: [it] }.combine(ch_intervals_sthresh).map{ it: it[0] }
        ch_intervals_gwa = GCTA_PERFORM_GWA.out.gwa.map{ it: [it] }.combine(ch_intervals_sthresh).map{ it: it[0] }

        R_GET_GCTA_INTERVALS(
            ch_intervals_params,
            ch_intervals_grm,
            ch_intervals_plink,
            ch_intervals_pheno,
            ch_intervals_gwa,
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
}


// Current bug that publish doesn't work without an output closure
output {
    "." {
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
    Strain Name and List File               = ${params.strainfile}
    VCF                                     = ${params.vcf}
    Number of simulated QTLs                = ${nqtl_file}
    Phenotype Heritability File             = ${h2_file}
    Number of simulation replicates         = ${params.reps}
    MAF Threshold File                      = ${maf_file}
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

    // println summary

}
