process PLINK_UPDATE_BY_H2 {

    label 'plink_update_by_h2'
    tag "${nqtl} ${rep} ${h2} ${effect} ${group}_${maf}"

    input:
    // filter_id (CV region label) rides as a trailing val so it reaches the GWA DB write.
    // pool_hash arrives from GCTA_SIMULATE_PHENOTYPES.out.params but is only needed on the
    // trait branch, so it is accepted here and dropped from the output.
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(species), val(filter_id), val(pool_hash)
    tuple val(group1), val(maf1), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests)
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par")

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed"),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim"),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam"),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map"),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex"),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped"),
          path("TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log"),
          path(gm), path(n_indep_tests),
          path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"),
          path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"),
          val(species), val(filter_id),
          emit: out
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

    plink --bfile TO_SIMS \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${maf} \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group} \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --pheno ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$( plink --version |& head -n 1 | cut -f 2 )
    END_VERSIONS
    """

    stub:
    """
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bed
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.fam
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.map
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.nosex
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.ped
    touch TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: stub
    END_VERSIONS
    """
}
