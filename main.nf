nextflow.enable.dsl=2

// Needed to publish results
nextflow.preview.output = true

// import the subworkflows
include { LOCAL_GET_CONTIG_INFO           } from './modules/local/get_contig_info/main'
include { LOCAL_COMPILE_EIGENS            } from './modules/local/compile_eigens/main'
include { BCFTOOLS_EXTRACT_STRAINS        } from './modules/bcftools/extract_strains/main'
include { BCFTOOLS_CREATE_GENOTYPE_MATRIX } from './modules/bcftools/create_genotype_matrix/main'
include { PLINK_RECODE_VCF                } from './modules/plink/recode_vcf/main'
include { PLINK_UPDATE_BY_H2              } from './modules/plink/update_by_h2/main'
include { R_FIND_GENOTYPE_MATRIX_EIGEN    } from './modules/r/find_genotype_matrix_eigen/main'
include { R_SIMULATE_EFFECTS_LOCAL        } from './modules/r/simulate_effects_local/main'
include { R_SIMULATE_EFFECTS_GLOBAL       } from './modules/r/simulate_effects_global/main'
include { R_GET_GCTA_INTERVALS           } from './modules/r/get_gcta_intervals/main'
include { R_ASSESS_SIMS                   } from './modules/r/assess_sims/main'
include { GCTA_SIMULATE_PHENOTYPES        } from './modules/gcta/simulate_phenotypes/main'
include { GCTA_MAKE_GRM                   } from './modules/gcta/make_grm/main'
include { GCTA_PERFORM_GWA                } from './modules/gcta/perform_gwa/main'

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
    log.info "--rep             Integer            The number of replicates to simulate per number of QTL and heritability (Default: 2)"
    log.info "--maf             File               A CSV file where each line is a minor allele frequency threshold to test for simulations (Default: data/simulate_maf.csv)"
    log.info "--effect          File               A CSV file where each line is an effect size range (e.g. 0.2-0.3) to test for simulations (Default: data/simulate_effect_sizes.csv)"
    log.info "--qtlloc          File               A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])"
    log.info "--sthresh         String             Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
    log.info "--group_qtl       Integer            If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
    log.info "--ci_size         Integer            Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
    log.info "--sparse_cut      Decimal            Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
    log.info "--mito_name       Strain             Name of mitochondrial chromosome"
    log.info "--simulate_qtlloc Boolean            Whether to simulate QTLs in specific genomic regions (Default: false)"
    log.info "-output-dir       String             Name of folder that will contain the results (Default: Simulations_{date})"
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
    log.info "Significance Threshold                  = ${params.sthresh}"
    log.info "Window for combining QTLs               = ${params.group_qtl}"
    log.info "Number of SNVs to define QTL CI         = ${params.ci_size}"
    log.info "Relatedness cutoff                      = ${params.sparse_cut}"
    log.info "Mitochondrial chromosome name           = ${mito_name}"
    log.info "Simulate QTLs in specific regions       = ${simulate_qtlloc}"
    log.info "Output directory                        = ${workflow.outputDir}"
    log.info ""
}



    // Created needed channels
    ch_vcf = Channel.fromPath(params.vcf).map{ it: [[id: "vcf"], it, "${it}.tbi"] }
    ch_strain_sets = Channel.fromPath(strainfile).splitCsv(sep: " ").map{ it: [[id: it[0]], it[1]] }
    ch_mafs = Channel.fromPath(maf_file).splitCsv().first()

    // Get contig data from VCF file
    LOCAL_GET_CONTIG_INFO( ch_vcf )
    ch_mito_num = LOCAL_GET_CONTIG_INFO.out.mapping.splitCsv(sep:"\t")
        .filter{ row -> row[0] == mito_name }
        .map{ row -> row[1] }
        .first()

    // Extract desired strain sets
    BCFTOOLS_EXTRACT_STRAINS( ch_vcf,
                              ch_strain_sets,
                              LOCAL_GET_CONTIG_INFO.out.mapping )
    ch_versions = ch_versions.mix(BCFTOOLS_EXTRACT_STRAINS.out.versions)

    // Recode the VCF file and create plink formatted files
    PLINK_RECODE_VCF( BCFTOOLS_EXTRACT_STRAINS.out.vcf,
                      ch_mito_num,
                      ch_mafs )
    ch_versions = ch_versions.mix(PLINK_RECODE_VCF.out.versions)

    // Create plaintext genotype matrix
    BCFTOOLS_CREATE_GENOTYPE_MATRIX( PLINK_RECODE_VCF.out.vcf,
                                     PLINK_RECODE_VCF.out.markers )
    ch_versions = ch_versions.mix(BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.versions)

    // Find eigen values for genotype matrix
    ch_chrom_nums = LOCAL_GET_CONTIG_INFO.out.mapping
        .splitCsv(sep: "\t")
        .filter{ it: it[0] != mito_name }
        .map{ it: it[1] }
        .toSortedList()

    ch_chrom_nums.view()

    R_FIND_GENOTYPE_MATRIX_EIGEN( BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix,
                                  Channel.fromPath("${workflow.projectDir}/bin/Get_GenoMatrix_Eigen.R").first(),
                                  ch_chrom_nums )
    ch_versions = ch_versions.mix(R_FIND_GENOTYPE_MATRIX_EIGEN.out.versions)

    // Concatenate eigen files with plink and genotype matrix sets
    ch_eigens = R_FIND_GENOTYPE_MATRIX_EIGEN.out.eigen
        .groupTuple(by: [0, 1])
        .map{ it: [it[0], it[1], it[3]] }

    LOCAL_COMPILE_EIGENS( ch_eigens )

    // Compile required files for simulations by strain group and MAF
    ch_plink_genomat_eigen = PLINK_RECODE_VCF.out.plink
        .join(BCFTOOLS_CREATE_GENOTYPE_MATRIX.out.matrix, by: [0, 1])
        .join(ch_eigens, by: [0, 1])
        .join(LOCAL_COMPILE_EIGENS.out.tests, by: [0, 1])
    
    // Simulate QTL or genome
    if (simulate_qtlloc){
        R_SIMULATE_EFFECTS_LOCAL( ch_plink_genomat_eigen,
                                  Channel.fromPath(params.qtlloc),
                                  Channel.fromPath("${workflow.projectDir}/bin/Create_Causal_QTLs.R").first(),
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
        ch_versions = ch_versions.mix(R_SIMULATE_EFFECTS_LOCAL.out.versions)
        ch_sim_phenos = R_SIMULATE_EFFECTS_LOCAL.out.causal
        ch_sim_plink = R_SIMULATE_EFFECTS_LOCAL.out.plink
    } else {
        R_SIMULATE_EFFECTS_GLOBAL( ch_plink_genomat_eigen,
                                   Channel.fromPath("${workflow.projectDir}/bin/Create_Causal_QTLs.R").first(),
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
        ch_versions = ch_versions.mix(R_SIMULATE_EFFECTS_GLOBAL.out.versions)
        ch_sim_phenos = R_SIMULATE_EFFECTS_GLOBAL.out.causal
        ch_sim_plink = R_SIMULATE_EFFECTS_GLOBAL.out.plink
    }

    // Adjust plink data by heritability
    GCTA_SIMULATE_PHENOTYPES( ch_sim_phenos,
                              ch_sim_plink,
                              Channel.fromPath("${h2_file}").splitCsv().first() )
    ch_versions = ch_versions.mix(GCTA_SIMULATE_PHENOTYPES.out.versions)

    // Update plink data by heritability
    PLINK_UPDATE_BY_H2( GCTA_SIMULATE_PHENOTYPES.out.params,
                        GCTA_SIMULATE_PHENOTYPES.out.plink,
                        GCTA_SIMULATE_PHENOTYPES.out.pheno )
    ch_versions = ch_versions.mix(PLINK_UPDATE_BY_H2.out.versions)

    // Create genetic relatedness matrix
    ch_mode = Channel.of( ["inbred", "fastGWA"], ["loco", "mlma"] )
    ch_grm_params = PLINK_UPDATE_BY_H2.out.params.combine(ch_mode)
    ch_grm_plink = PLINK_UPDATE_BY_H2.out.plink.map{ it: [it] }.combine(ch_mode).map{ it: it[0] }
    ch_grm_pheno = PLINK_UPDATE_BY_H2.out.pheno.map{ it: [it] }.combine(ch_mode).map{ it: it[0] }

    ch_grm_params.view()
    ch_grm_plink.view()
    ch_grm_pheno.view()

    GCTA_MAKE_GRM( ch_grm_params,
                   ch_grm_plink,
                   ch_grm_pheno )
    ch_versions = ch_versions.mix(GCTA_MAKE_GRM.out.versions)

    // Simulate GWA
    ch_type = Channel.of( "pca", "nopca" )
    ch_gwa_params = GCTA_MAKE_GRM.out.params.combine(ch_type)
    ch_gwa_grm = GCTA_MAKE_GRM.out.plink.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    ch_gwa_plink = GCTA_MAKE_GRM.out.plink.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    ch_gwa_pheno = GCTA_MAKE_GRM.out.pheno.map{ it: [it] }.combine(ch_type).map{ it: it[0] }



    GCTA_PERFORM_GWA( ch_gwa_params,
                      ch_gwa_grm,
                      ch_gwa_plink,
                      ch_gwa_pheno,
                      params.sparse_cut )
    ch_versions = ch_versions.mix(GCTA_PERFORM_GWA.out.versions)

    // Find GCTA intervals
    R_GET_GCTA_INTERVALS( GCTA_PERFORM_GWA.out.params,
                           GCTA_PERFORM_GWA.out.grm,
                           GCTA_PERFORM_GWA.out.plink,
                           GCTA_PERFORM_GWA.out.pheno,
                           GCTA_PERFORM_GWA.out.gwa,
                           Channel.fromPath("${workflow.projectDir}/bin/Get_GCTA_Intervals.R").first(),
                           params.sthresh,
                           params.group_qtl,
                           params.ci_size )
    ch_versions = ch_versions.mix(R_GET_GCTA_INTERVALS.out.versions)

    // Compile results
    R_ASSESS_SIMS( R_GET_GCTA_INTERVALS.out.params,
                   R_GET_GCTA_INTERVALS.out.grm,
                   R_GET_GCTA_INTERVALS.out.plink,
                   R_GET_GCTA_INTERVALS.out.pheno,
                   R_GET_GCTA_INTERVALS.out.interval,
                   Channel.fromPath("${workflow.projectDir}/bin/Assess_Sims.R").first() )
    ch_versions = ch_versions.mix(R_ASSESS_SIMS.out.versions)

    // Split results by algorithm and compile into summary file
    R_ASSESS_SIMS.out.assessment.branch{ v ->
        inbred: v.contains("inbred_nopca")
        inbred_pca: v.contains("inbred_pca")
        loco: v.contains("loco_nopca")
        loco_pca: v.contains("loco_pca")
        }.set{ result }

    publish:
    result.inbred     >> "."
    result.inbred_pca >> "."
    result.loco       >> "."
    result.loco_pca   >> "."
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
    Significance Threshold                  = ${params.sthresh}
    Threshold for grouping QTL              = ${params.group_qtl}
    Number of SNVs to define CI             = ${params.ci_size}
    Relatedness Matrix Cutoff               = ${params.sparse_cut}
    Mitochondrial chromosome name           = ${mito_name}
    Simulate QTLs in specific regions       = ${simulate_qtlloc}
    Result Directory                        = ${workflow.outputDir}
    """

    // println summary

}
