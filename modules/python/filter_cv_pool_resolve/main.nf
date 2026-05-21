process FILTER_CV_POOL_RESOLVE {

    label 'python_filter_cv_pool_resolve'
    tag "${group}_${ms_maf}_${filter_id}"

    input:
    tuple val(group), val(ms_maf),
          path("CV_TO_SIMS.bed"), path("CV_TO_SIMS.bim"), path("CV_TO_SIMS.fam"),
          val(species), val(filter_id)
    path nqtl_file
    val reps
    val rep_start
    path domain_table        // per-species recombination-domain reference (Change 3)
    path filter_script

    output:
    tuple val(group), val(ms_maf), val(filter_id),
          path("CV_TO_SIMS.bed"), path("CV_TO_SIMS.bim"), path("CV_TO_SIMS.fam"),
          path("extract.list"),
          path("resolved_intervals.json"),     // sidecar (Change 3)
          path("pool_hash.txt"),               // sha256(resolved marker list) -- single-source seed
        emit: extract_list
    path "filter_metrics.tsv", emit: metrics   // one row per nqtl
    path "versions.yml",       emit: versions

    script:
    """
    python ${filter_script} \\
        --bim CV_TO_SIMS.bim \\
        --filter_id ${filter_id} \\
        --species ${species} \\
        --domain_file ${domain_table} \\
        --nqtl_file ${nqtl_file} \\
        --reps ${reps} \\
        --rep_start ${rep_start} \\
        --group ${group} \\
        --ms_maf ${ms_maf} \\
        --extract_out extract.list \\
        --resolved_out resolved_intervals.json \\
        --pool_hash_out pool_hash.txt \\
        --metrics_out filter_metrics.tsv

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$( python --version 2>&1 | cut -d' ' -f2 )
            numpy: \$( python -c "import numpy; print(numpy.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    cp CV_TO_SIMS.bim extract.list
    echo '[]' > resolved_intervals.json
    echo stub > pool_hash.txt
    echo -e "group\\tms_maf\\tfilter_id\\tspecies\\tnqtl\\tpool_size\\tn_combinations\\trep_start\\trequested_reps\\teffective_reps\\tstatus" > filter_metrics.tsv
    while read nq; do
        [ -z "\$nq" ] && continue
        echo -e "${group}\\t${ms_maf}\\t${filter_id}\\t${species}\\t\${nq}\\t0\\t0\\t${rep_start}\\t${reps}\\t0\\tstarved" >> filter_metrics.tsv
    done < ${nqtl_file}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: stub
    END_VERSIONS
    """
}
