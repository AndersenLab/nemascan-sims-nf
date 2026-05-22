process GCTA_SIMULATE_PHENOTYPES {

    label 'gcta_simulate_phenotypes'
    tag "${species}_${group}_${maf}_${nqtl}_${effect}_${h2}_${rep}"

    input:
    // filter_id (CV region label) + pool_hash (sha of the resolved CV pool) ride as trailing
    // vals; the trait write uses them to record the region and derive the causal-set id in R.
    // CV pool is bed/bim/fam only — the upstream region-filter extraction emits just those.
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), path(causal_variants), val(species), val(filter_id), val(pool_hash)
    tuple val(group), val(maf),
          path("CV_TO_SIMS.bed"),  path("CV_TO_SIMS.bim"),  path("CV_TO_SIMS.fam")
    tuple val(group), val(maf),
          path("TO_SIMS.bed"),     path("TO_SIMS.bim"),      path("TO_SIMS.fam"),
          path("TO_SIMS.map"),     path("TO_SIMS.nosex"),    path("TO_SIMS.ped"),    path("TO_SIMS.log"),
          path(gm), path(n_indep_tests)

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(species), val(filter_id), val(pool_hash), emit: params
    tuple val(group), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests), emit: plink
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"), emit: pheno
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export NF_TRAP_FAILURES_DIR="${workflow.outputDir}/.failures"
    export NF_TRAP_CELL_KEY="${species}__${group}__${maf}__${nqtl}__${effect}__${rep}__${h2}__${filter_id}__NA__NA"
    export NF_TRAP_PAYLOAD='{"session":"${workflow.sessionId}","task_hash":"${task.hash}","attempt":${task.attempt},"max_retries":${task.maxRetries},"species":"${species}","group":"${group}","maf":${maf},"nqtl":"${nqtl}","effect":"${effect}","h2":${h2},"rep":${rep},"filter_id":"${filter_id}","mode":"NA","type":"NA"}'
    source ${projectDir}/bin/failure_trap.sh

    gcta64 --bfile CV_TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${causal_variants} \\
         --simu-hsq ${h2} \\
         --simu-rep ${rep} \\
         --autosome-num 6 \\
         --thread-num 1 \\
         --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims  # pinned: BLAS reduction order must be deterministic

    awk -v col=\$((${rep} + 2)) '{print \$1, \$2, \$col}' \\
        ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen \\
        > ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen.tmp && \\
        mv ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen.tmp \\
           ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """
    
    stub:
    """
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par
    touch CV_TO_SIMS.bed
    touch CV_TO_SIMS.bim
    touch CV_TO_SIMS.fam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: stub
    END_VERSIONS
    """
}
