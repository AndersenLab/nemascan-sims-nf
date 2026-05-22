process DB_MIGRATION_WRITE_GWA_TO_DB {

    label 'db_migration_write_gwa_to_db'
    tag "${species}_${group}_${maf}_${nqtl}_${effect}_${h2}_${rep}_${mode}_${type}"

    input:
    // filter_id (CV region label) rides as a trailing val; write_gwa_to_db.R re-derives the
    // region-aware trait_id → mapping_id in R from it.
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2),
          val(mode), val(type),
          val(cv_maf_effective), val(cv_ld), val(species), val(filter_id)
    path gwa_file
    val base_dir

    output:
    val true, emit: done
    path "versions.yml", emit: versions

    script:
    """
    export NF_TRAP_FAILURES_DIR="${workflow.outputDir}/.failures"
    export NF_TRAP_CELL_KEY="${species}__${group}__${maf}__${nqtl}__${effect}__${rep}__${h2}__${filter_id}__${mode}__${type}"
    export NF_TRAP_PAYLOAD='{"session":"${workflow.sessionId}","task_hash":"${task.hash}","attempt":${task.attempt},"max_retries":${task.maxRetries},"species":"${species}","group":"${group}","maf":${maf},"nqtl":"${nqtl}","effect":"${effect}","h2":${h2},"rep":${rep},"filter_id":"${filter_id}","mode":"${mode}","type":"${type}"}'
    source ${projectDir}/bin/failure_trap.sh

    export R_SOURCE_DIR="${projectDir}/R"
    write_gwa_to_db.R \
        --group ${group} --maf ${maf} --nqtl ${nqtl} --effect ${effect} \
        --rep ${rep} --h2 ${h2} --mode ${mode} --type ${type} \
        --cv_maf_effective ${cv_maf_effective} --cv_ld ${cv_ld} \
        --cv_region_filter ${filter_id} \
        --gwa_file ${gwa_file} --base_dir ${base_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    echo "STUB: write_gwa_to_db ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_${mode}_${type} cv_maf=${cv_maf_effective} cv_ld=${cv_ld}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
