process GCTA_SIMULATE_PHENOTYPES {

    label 'gcta_simulate_phenotypes'
    tag "${nqtl} ${rep} ${h2} ${effect} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), path(causal_variants)
    tuple val(group), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests)
    each h2

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), emit: params
    tuple val(group), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests), emit: plink
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"), emit: pheno
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gcta64 --bfile TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${causal_variants} \\
         --simu-hsq ${h2} \\
         --simu-rep 1 \\
         --autosome-num 6 \\
         --thread-num ${task.cpus} \\
         --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """

    stub:
    """
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """
}
