process DB_MIGRATION_WRITE_MARKER_SET {

    label 'db_migration_write_marker_set'
    tag "${group}_${maf}"

    input:
    tuple val(group), val(maf), path(bim_file), path(n_indep_tests)
    val base_dir

    output:
    val true, emit: done

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    write_marker_set.R --group ${group} --maf ${maf} --bim ${bim_file} --eigen ${n_indep_tests} --base_dir ${base_dir}
    """

    stub:
    """
    echo "STUB: write_marker_set ${group}_${maf}"
    """
}
