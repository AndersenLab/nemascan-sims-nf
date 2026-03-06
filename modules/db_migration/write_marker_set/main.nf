process DB_MIGRATION_WRITE_MARKER_SET {

    label 'db_migration_write_marker_set'
    tag "${group}_${maf}"

    input:
    tuple val(group), val(maf), path(bim_file), path(n_indep_tests),
          val(species), val(vcf_release_id), val(ms_ld), val(strains), val(strainfile_path)
    val base_dir

    output:
    val true, emit: done
    path "versions.yml", emit: versions

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    write_marker_set.R --group ${group} --maf ${maf} --bim ${bim_file} --eigen ${n_indep_tests} \
        --base_dir ${base_dir} \
        --species ${species} \
        --vcf_release_id ${vcf_release_id} \
        --ms_ld ${ms_ld}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
        --ms_ld ${ms_ld} \
        --strainfile_path "${strainfile_path}" \
        --strains "${strains}"
    """

    stub:
    """
    echo "STUB: write_marker_set ${group}_${maf} species=${species} vcf_release_id=${vcf_release_id} ms_ld=${ms_ld}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
