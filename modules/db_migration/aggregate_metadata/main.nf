process DB_MIGRATION_AGGREGATE_METADATA {

    label 'db_migration_aggregate_metadata'
    tag "aggregate"

    input:
    val ready
    val base_dir

    output:
    path "aggregation_summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    aggregate_metadata.R ${base_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    echo "STUB: aggregate_metadata" > aggregation_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
