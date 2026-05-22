process VALIDATE_REPLICATION_COMPLETE {
    label 'r_validate_replication'
    tag 'validation'

    input:
    val trigger
    val db_root
    val expected_reps
    path filter_pool_metrics

    output:
    val true, emit: done

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    validate_replication_complete.R \
        --db_root ${db_root} \
        --expected_reps ${expected_reps} \
        --filter_pool_metrics ${filter_pool_metrics} \
        --output_dir "${workflow.outputDir}"
    """

    stub:
    """
    echo "STUB: validate_replication_complete (db_root=${db_root}, reps=${expected_reps}, metrics=${filter_pool_metrics})"
    """
}
