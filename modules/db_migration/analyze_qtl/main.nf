process DB_MIGRATION_ANALYZE_QTL {

    label 'db_migration_analyze_qtl'
    tag "${threshold} ${nqtl} ${rep} ${h2} ${effect} ${mode} ${type} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), val(threshold)
    tuple path(pheno_file), path(par_file)
    val base_dir
    val ci_size
    val snp_grouping
    val alpha

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), val(threshold), emit: params
    tuple path(pheno_file), path(par_file), emit: pheno
    path "*_qtl_regions.tsv", emit: regions
    val true, emit: done
    path "versions.yml", emit: versions

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    analyze_qtl.R \
        --group ${group} --maf ${maf} --nqtl ${nqtl} --effect ${effect} \
        --rep ${rep} --h2 ${h2} --mode ${mode} --type ${type} \
        --threshold ${threshold} \
        --base_dir ${base_dir} --ci_size ${ci_size} --snp_grouping ${snp_grouping} \
        --alpha ${alpha}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    touch "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_${mode}_${type}_${threshold}_qtl_regions.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
