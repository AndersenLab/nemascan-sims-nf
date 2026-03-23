process DB_MIGRATION_WRITE_GWA_TO_DB {

    label 'db_migration_write_gwa_to_db'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${type} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2),
          val(mode), val(type),
          val(cv_maf_effective), val(cv_ld)
    path gwa_file
    val base_dir

    output:
    val true, emit: done
    path "versions.yml", emit: versions

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    write_gwa_to_db.R \
        --group ${group} --maf ${maf} --nqtl ${nqtl} --effect ${effect} \
        --rep ${rep} --h2 ${h2} --mode ${mode} --type ${type} \
        --cv_maf_effective ${cv_maf_effective} --cv_ld ${cv_ld} \
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
