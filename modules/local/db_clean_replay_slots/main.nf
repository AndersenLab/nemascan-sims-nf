process DB_CLEAN_REPLAY_SLOTS {
    label 'db_migration_cleanup'
    tag 'replay-cleanup'

    errorStrategy { task.attempt <= task.maxRetries ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path replay_tsv
    val  db_root
    val  cv_maf
    val  cv_ld

    output:
    val 'done', emit: done

    script:
    def cv_maf_arg = cv_maf != null ? "--cv_maf ${cv_maf}" : ""
    """
    export R_SOURCE_DIR="${projectDir}/R"
    db_clean_replay_slots.R \\
        --replay_tsv ${replay_tsv} \\
        --db_root ${db_root} \\
        ${cv_maf_arg} \\
        --cv_ld ${cv_ld}
    """

    stub:
    """
    echo "STUB: db_clean_replay_slots ${replay_tsv} → ${db_root} cv_maf=${cv_maf} cv_ld=${cv_ld}"
    """
}
