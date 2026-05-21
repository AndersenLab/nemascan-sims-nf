process FILTER_CV_POOL_EXTRACT {

    label 'plink_recode_vcf'   // reuse existing PLINK label + container
    tag "${group}_${ms_maf}_${filter_id}"

    input:
    tuple val(group), val(ms_maf), val(filter_id),
          path("CV_TO_SIMS.bed"), path("CV_TO_SIMS.bim"), path("CV_TO_SIMS.fam"),
          path("extract.list"),
          path("resolved_intervals.json"),
          path("pool_hash.txt")

    output:
    tuple val(group), val(ms_maf), val(filter_id),
          path("CV_FILTERED.bed"), path("CV_FILTERED.bim"), path("CV_FILTERED.fam"),
          path("resolved_intervals.json"),        // sidecar passthrough
          path("pool_hash.txt"),                  // passthrough (single-source seed)
        emit: plink
    path "versions.yml", emit: versions

    script:
    """
    plink --bfile CV_TO_SIMS \\
        --extract extract.list \\
        --make-bed \\
        --out CV_FILTERED \\
        --allow-extra-chr

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            plink: \$( plink --version | head -1 | awk '{print \$2}' )
    END_VERSIONS
    """

    stub:
    """
    for ext in bed bim fam; do
        touch CV_FILTERED.\$ext
    done

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            plink: stub
    END_VERSIONS
    """
}
