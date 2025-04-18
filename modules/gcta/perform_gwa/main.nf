process GCTA_PERFORM_GWA {

    label 'gcta_perform_gwa'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${type} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type) 
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id")
    tuple path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests), path(causal_variants)
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par")
    val sparse_cut

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type),                                     emit: params
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id"), emit: grm
    tuple path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIM_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests), path(causal_variants), emit: plink
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"),     emit: pheno
    path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_${mode}_${type}.${suffix}",                                      emit: gwa
    path "versions.yml",                                                                                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [[ ${mode} == "inbred" ]]; then
        GRM_OPTION="--grm-sparse"
        COMMAND="--fastGWA-lmm-exact"
    else
        GRM_OPTION="--grm"
        COMMAND="--mlma-loco"
    fi

    if [[ ${type} == "pca" ]]; then
        gcta64 --grm TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --make-bK-sparse ${sparse_cut} \\
            --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sparse_grm_${mode} \\
            --thread-num ${task.cpus}

        gcta64 --grm TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --pca 1 \\
            --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sparse_grm_${mode} \\
            --thread-num ${task.cpus}

        COVAR="--qcovar ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sparse_grm_${mode}.eigenvec"
    else
        COVAR=""
    fi

    gcta64 \${COMMAND} \\
        --bfile TO_SIMS_${nqtl}_${rep}_${maf}_${effect}_${group} \\
        \${GRM_OPTION} ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sparse_grm_${mode} \\
        \${COVAR} \\
        --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type} \\
        --pheno ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen \\
        --maf ${maf} \\
        --thread-num ${task.cpus}

    if [[ ${mode} == "loco" ]]; then
        mv "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_${mode}_${type}.loco.mlma" \\
           "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_${mode}_${type}.mlma"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """

    stub:
    """
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """
}