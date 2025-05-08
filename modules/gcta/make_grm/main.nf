process GCTA_MAKE_GRM {

    label 'gcta_make_grm'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${group}_${maf}"

    input:

    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix)
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests)
    tuple path(pheno_in_tmp), path(pheno_in_par)

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), emit: params
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id"), emit: grm
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests), emit: plink
    tuple path(pheno_in_tmp.name), path("check_vp.hsq"), path(pheno_in_par), emit: pheno_hsq_and_par
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [[ ${mode} == "inbred" ]]; then
        GRM_OPTION="--make-grm-inbred"
    else
        GRM_OPTION="--make-grm"
    fi

    gcta64 --bfile TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group} \\
            --autosome --maf ${maf} \${GRM_OPTION} \\
            --out TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --thread-num ${task.cpus}

    gcta64 --grm TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --pheno ${pheno_in_tmp} \\
            --reml --out check_vp \\
            --thread-num ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """

    stub:
    """
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id
    touch "${pheno_in_tmp.name}"
    touch "check_vp.hsq"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """
}
