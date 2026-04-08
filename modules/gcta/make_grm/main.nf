process GCTA_MAKE_GRM {

    label 'gcta_make_grm'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${group}_${maf}"

    input:

    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix)
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests)
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par")

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), emit: params
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.bin"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N.bin"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id"), emit: grm
    tuple path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"), path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"), path(gm), path(n_indep_tests), emit: plink
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"), emit: pheno
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

    awk '{print \$2}' TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim > plink_snplist.txt

    gcta64 --bfile TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group} \\
            --autosome --extract plink_snplist.txt \${GRM_OPTION} \\
            --out TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode} \\
            --thread-num 1  # pinned: BLAS reduction order must be deterministic

    # Iterative REML-scale loop: run REML, check Vp, scale x1000 per round
    # until full-GRM Vp >= 1e-4 (max 4 rounds).
    VP_THRESHOLD="0.0001"
    VP_MAX_ROUNDS=4
    PHEN_INPUT="${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"
    GRM_PREFIX="TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}"

    cp -L "\${PHEN_INPUT}" working_pheno.txt   # copy from symlink

    for vp_round in \$(seq 1 \${VP_MAX_ROUNDS}); do
        if ! gcta64 --grm \${GRM_PREFIX} \\
                    --pheno working_pheno.txt \\
                    --reml --out check_vp \\
                    --thread-num 1; then
            echo "REML failed (round \${vp_round}); falling back to Vp=1.0"
            printf 'Source\\tVariance\\tSE\\nVp\\t1.0\\tNA\\n' > check_vp.hsq
            break
        fi

        VP_VALUE=\$(awk -F'\\t' '\$1 == "Vp" {print \$2}' check_vp.hsq)
        VP_OK=\$(awk -v vp="\${VP_VALUE}" -v thresh="\${VP_THRESHOLD}" \\
                 'BEGIN { print (vp+0 >= thresh+0) ? "1" : "0" }')

        if [ "\${VP_OK}" -eq 1 ]; then
            echo "Vp=\${VP_VALUE} >= \${VP_THRESHOLD} at round \${vp_round}; done"
            break
        fi

        if [ "\${vp_round}" -eq "\${VP_MAX_ROUNDS}" ]; then
            echo "WARNING: Vp=\${VP_VALUE} still below \${VP_THRESHOLD} after \${VP_MAX_ROUNDS} rounds"
            break
        fi

        echo "Vp=\${VP_VALUE} < \${VP_THRESHOLD} at round \${vp_round}; scaling x1000"
        awk '{printf "%s %s %.10g\\n", \$1, \$2, \$3 * 1000}' working_pheno.txt > scaled_pheno.tmp
        mv scaled_pheno.tmp working_pheno.txt
    done

    cp working_pheno.txt "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: \$( gcta64 --version |& grep version |& cut -f 3 )
    END_VERSIONS
    """

    stub:
    """
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.bin
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.N.bin
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}.grm.id
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno
    touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GCTA: stub
    END_VERSIONS
    """
}
