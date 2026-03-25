process GCTA_PERFORM_GWA {

    label 'gcta_perform_gwa'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${type} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type)
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.bin"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N.bin"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id")
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests)
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par")
    val sparse_cut

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), emit: params
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.bin"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N.bin"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id"), emit: grm
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests), emit: plink
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"), emit: pheno
    path "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type}.${suffix}", emit: gwa
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [[ ${mode} == "inbred" ]]; then
        COMMAND='--fastGWA-mlm-exact'
        GWA_THREADS=1  # pinned: BLAS reduction order must be deterministic

        # Inbred: create sparse GRM approximation (supported by fastGWA-mlm-exact)
        gcta64 --grm TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --make-bK-sparse ${sparse_cut} \\
            --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sparse_grm_${mode} \\
            --thread-num 1

        GRM_OPTION='--grm-sparse'
        GRM_PREFIX=${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sparse_grm_${mode}
    else
        COMMAND="--mlma-loco"
        GWA_THREADS=${task.cpus}  # mlma-loco: per-chromosome refits are independent

        # LOCO: must use the full (non-sparse) GRM — --mlma-loco decomposes per-chromosome
        # covariance from the full GRM; a sparse approximation causes REML divergence.
        GRM_OPTION="--grm"
        GRM_PREFIX=TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}
    fi

    if [[ ${type} == "pca" ]]; then

        gcta64 --grm TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --pca 1 \\
            --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_pca \\
            --thread-num 1  # pinned: BLAS reduction order must be deterministic

        COVAR="--qcovar ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_pca.eigenvec"
    else
        COVAR=""
    fi

    awk '{print \$2}' TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim > plink_snplist.txt

    gcta64 \${COMMAND} \\
        --bfile TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group} \\
        \${GRM_OPTION} \${GRM_PREFIX} \\
        \${COVAR} \\
        --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type} \\
        --pheno ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno \\
        --extract plink_snplist.txt \\
        --thread-num \${GWA_THREADS}

    if [[ ${mode} == "loco" ]]; then
        mv "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type}.loco.mlma" \\
           "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type}.mlma"
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
        GCTA: stub
    END_VERSIONS
    """
}
