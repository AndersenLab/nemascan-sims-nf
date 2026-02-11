process DB_MIGRATION_AGGREGATE_METADATA {

    label 'db_migration_aggregate_metadata'
    tag "aggregate"

    input:
    val ready
    val base_dir

    output:
    path "aggregation_summary.txt", emit: summary

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    aggregate_metadata.R ${base_dir}
    """

    stub:
    """
    echo "STUB: aggregate_metadata" > aggregation_summary.txt
    """
}
