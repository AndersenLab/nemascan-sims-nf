process BCFTOOLS_EXTRACT_STRAINS {

    label "bcftools_extract_strains"
    tag "${meta1}"

    input:
    tuple val(meta), path(vcf_renamed), path(vcf_renamed_index)
    tuple val(meta1), val(strains)

    output:
    tuple val(meta1), path("${meta1}.vcf.gz"), path("${meta1}.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools view --threads ${task.cpus} -s ${strains} -O u ${vcf_renamed} | \\
    bcftools filter -i N_MISSING=0 -O z -o ${meta1}.vcf.gz
    bcftools index --tbi ${meta1}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta1}.vcf.gz
    touch ${meta1}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}