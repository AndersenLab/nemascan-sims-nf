process R_ASSESS_SIMS {

    label 'r_assess_sims'
    tag "${group}_${maf}_${nqtl}_${effect}_${h2}_${rep}_${mode}_${type}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), val(threshold)
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id")
    tuple path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests)
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par")
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_processed_LMM-EXACT_${mode}_${type}_mapping.tsv"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_LMM-EXACT_${mode}_${type}_qtl_region.tsv")
    path assess_sims
    val alpha
    val ci_size
    val snp_grouping

    output:
    path "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_${mode}_${type}_${threshold}_mapping.tsv", emit: assessment
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export NF_TRAP_SESSION_ID="${workflow.sessionId}"
    export NF_TRAP_FAILURES_DIR="${workflow.outputDir}/.failures"
    export NF_TRAP_TASK_HASH="${task.hash}"
    export NF_TRAP_ATTEMPT=${task.attempt}
    export NF_TRAP_MAX_RETRIES=${task.maxRetries}
    export GROUP="${group}"
    export MAF=${maf}
    export NQTL="${nqtl}"
    export EFFECT="${effect}"
    export H2=${h2}
    export REP=${rep}
    export MODE="${mode}"
    export TYPE="${type}"
    source ${projectDir}/templates/failure_trap.sh

    Rscript --vanilla ${assess_sims} ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_processed_LMM-EXACT_${mode}_${type}_mapping.tsv \\
        ${gm} ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno \\
        ${nqtl} ${rep} ${h2} ${maf} ${effect} ${group} ${mode}_${type}_${threshold} \\
        ${alpha} ${ci_size} ${snp_grouping}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    touch "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_${mode}_${type}_${threshold}_mapping.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
