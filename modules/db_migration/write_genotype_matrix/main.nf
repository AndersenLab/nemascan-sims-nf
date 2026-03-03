process DB_MIGRATION_WRITE_GENOTYPE_MATRIX {
    label 'db_migration_write_genotype_matrix'
    tag "${group}_${maf}"

    input:
    tuple val(group), val(maf), path(genotype_matrix),
          val(species), val(vcf_release_id), val(ms_ld)
    val base_dir

    output:
    val true, emit: done

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    write_genotype_matrix.R --group ${group} --maf ${maf} \
        --genotype_matrix ${genotype_matrix} --base_dir ${base_dir} \
        --species ${species} \
        --vcf_release_id ${vcf_release_id} \
        --ms_ld ${ms_ld}
    """

    stub:
    """
    echo "STUB: write_genotype_matrix ${group}_${maf} species=${species} vcf_release_id=${vcf_release_id} ms_ld=${ms_ld}"
    """
}
