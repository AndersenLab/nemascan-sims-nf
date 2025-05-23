process BCFTOOLS_RENAME_CHROMS {

    label "bcftools_rename_chroms"
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path "chr_mapping.tsv"

    output:
    tuple val(meta), path("${meta.id}_renamed.vcf.gz"), path("${meta.id}_renamed.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools view --threads ${task.cpus} -O u ${vcf} | \\
    bcftools annotate --rename-chrs chr_mapping.tsv -O z -o ${meta.id}_renamed.vcf.gz
    bcftools index --tbi ${meta.id}_renamed.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_renamed.vcf.gz
    touch ${meta.id}_renamed.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}