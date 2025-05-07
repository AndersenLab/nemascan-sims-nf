process BCFTOOLS_EXTRACT_STRAINS {

    label "bcftools_extract_strains"
    tag "${meta1.id}"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(meta1), val(strains)
    path "chr_mapping.tsv"

    output:
    tuple val(meta1), path("${meta1.id}.vcf.gz"), path("${meta1.id}.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools view --threads ${task.cpus} -s ${strains} -O u ${vcf} | \\
    bcftools annotate --rename-chrs chr_mapping.tsv - | \\
    bcftools filter -i N_MISSING=0 -O z -o ${meta1.id}.vcf.gz
    bcftools index --tbi ${meta1.id}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta1.id}.vcf.gz
    touch ${meta1.id}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}