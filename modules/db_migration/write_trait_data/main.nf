process DB_MIGRATION_WRITE_TRAIT_DATA {
    label 'db_migration_write_trait_data'
    tag "${group}_${nqtl}_${rep}_${h2}"

    input:
    // filter_id (CV region label) + pool_hash ride in the params tuple. write_trait_data.R
    // derives the region-filter hash and the causal_set_id in R from these.
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(filter_id), val(pool_hash)
    path pheno_file
    path par_file
    val base_dir
    path causal_geno_file
    val cv_maf_effective
    val cv_ld

    output:
    val true, emit: done
    path "versions.yml", emit: versions

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    write_trait_data.R --group ${group} --maf ${maf} --nqtl ${nqtl} \
        --effect ${effect} --rep ${rep} --h2 ${h2} \
        --pheno_file ${pheno_file} --par_file ${par_file} --base_dir ${base_dir} \
        --causal_geno_file ${causal_geno_file} \
        --cv_maf_effective ${cv_maf_effective} \
        --cv_ld ${cv_ld} \
        --cv_region_filter ${filter_id} \
        --pool_hash ${pool_hash}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    echo "STUB: write_trait_data ${group}_${nqtl}_${rep}_${h2}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
