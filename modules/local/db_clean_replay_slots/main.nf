process DB_CLEAN_REPLAY_SLOTS {
    label 'db_migration_cleanup'
    tag 'replay-cleanup'

    errorStrategy { task.attempt <= task.maxRetries ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path replay_tsv
    val  db_root

    output:
    val 'done', emit: done

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    db_clean_replay_slots.R \
        --replay_tsv ${replay_tsv} \
        --db_root ${db_root}
    """

    stub:
    """
    echo "STUB: db_clean_replay_slots ${replay_tsv} → ${db_root}"
    """
}
