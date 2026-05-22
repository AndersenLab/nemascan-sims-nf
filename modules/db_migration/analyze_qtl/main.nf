process DB_MIGRATION_ANALYZE_QTL {

    label 'db_migration_analyze_qtl'
    tag "${species}_${group}_${maf}_${nqtl}_${effect}_${h2}_${rep}_${mode}_${type}"

    input:
    // filter_id (CV region label) rides as a trailing val; forwarded so ASSESS_SIMS and
    // analyze_qtl.R can re-derive the region-aware trait_id in R.
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), val(threshold), val(cv_maf_effective), val(cv_ld), val(species), val(filter_id)
    tuple path(pheno_file), path(par_file)
    val base_dir
    val ci_size
    val snp_grouping
    val alpha

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), val(threshold), val(cv_maf_effective), val(cv_ld), val(species), val(filter_id), emit: params
    tuple path(pheno_file), path(par_file), emit: pheno
    path "*_qtl_regions.tsv", emit: regions
    val true, emit: done
    path "versions.yml", emit: versions

    script:
    """
    export NF_TRAP_FAILURES_DIR="${workflow.outputDir}/.failures"
    export NF_TRAP_CELL_KEY="${species}__${group}__${maf}__${nqtl}__${effect}__${rep}__${h2}__${filter_id}__${mode}__${type}"
    export NF_TRAP_PAYLOAD='{"session":"${workflow.sessionId}","task_hash":"${task.hash}","attempt":${task.attempt},"max_retries":${task.maxRetries},"species":"${species}","group":"${group}","maf":${maf},"nqtl":"${nqtl}","effect":"${effect}","h2":${h2},"rep":${rep},"filter_id":"${filter_id}","mode":"${mode}","type":"${type}"}'
    source ${projectDir}/bin/failure_trap.sh

    export R_SOURCE_DIR="${projectDir}/R"
    analyze_qtl.R \
        --group ${group} --maf ${maf} --nqtl ${nqtl} --effect ${effect} \
        --rep ${rep} --h2 ${h2} --mode ${mode} --type ${type} \
        --threshold ${threshold} \
        --cv_maf_effective ${cv_maf_effective} --cv_ld ${cv_ld} \
        --cv_region_filter ${filter_id} \
        --base_dir ${base_dir} --ci_size ${ci_size} --snp_grouping ${snp_grouping} \
        --alpha ${alpha}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    touch "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_${mode}_${type}_${threshold}_qtl_regions.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
