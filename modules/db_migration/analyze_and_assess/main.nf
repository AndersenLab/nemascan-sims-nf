process DB_MIGRATION_ANALYZE_AND_ASSESS {

    label 'db_migration_analyze_and_assess'
    tag "${threshold} ${nqtl} ${rep} ${h2} ${effect} ${mode} ${type} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), val(type), val(threshold)
    tuple path(pheno_file), path(par_file)
    val base_dir
    val ci_size
    val snp_grouping

    output:
    path "*_db_assessment.tsv", emit: assessment
    val true, emit: done

    script:
    """
    export R_SOURCE_DIR="${projectDir}/R"
    analyze_and_assess.R \
        --group ${group} --maf ${maf} --nqtl ${nqtl} --effect ${effect} \
        --rep ${rep} --h2 ${h2} --mode ${mode} --type ${type} \
        --threshold ${threshold} --par_file ${par_file} \
        --base_dir ${base_dir} --ci_size ${ci_size} --snp_grouping ${snp_grouping}
    """

    stub:
    """
    touch "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_${mode}_${type}_${threshold}_db_assessment.tsv"
    """
}
