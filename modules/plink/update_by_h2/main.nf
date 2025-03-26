process PLINK_UPDATE_BY_H2 {

    label 'plink_update_by_h2'
    tag "${nqtl} ${rep} ${h2} ${effect} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2)
    tuple(val(group), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"),
          path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests))
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par")

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2),                                                                        emit: params
    tuple(path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"),
          path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"),
          path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"),
          path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests), path(causal_variants)),               emit: plink
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"),     emit: pheno
    path "versions.yml",                                                                                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plink --bfile TO_SIMS \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${maf} \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out TO_SIMS_${nqtl}_${rep}_${maf}_${effect}_${group} \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --pheno ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$( plink --version |& head -n 1 | cut -f 2' )
    END_VERSIONS
    """

    stub:
    """
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped
    touch TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$( plink --version |& head -n 1 | cut -f 2' )
    END_VERSIONS
    """
}