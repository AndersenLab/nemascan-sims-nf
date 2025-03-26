process BCFTOOLS_CREATE_GENOTYPE_MATRIX {

    label "bcftools_create_genotype_matrix"
    tag "${group}_${maf}"

    input:
    tuple val(group), val(maf), path(vcf), path(vcf_index)
    tuple val(group), val(maf), path(markers)

    output:
    tuple val(group), val(maf), path("${group}_${maf}_Genotype_Matrix.tsv"), emit: matrix
    path "versions.yml",                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools query -l ${vcf} | \\
        sort > sorted_samples.txt

    bcftools view -v snps -S sorted_samples.txt -R ${markers} ${vcf} | \\
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' | \\
        sed 's/[[# 0-9]*]//g' | \\
        sed 's/:GT//g' | \\
        sed 's/0|0/-1/g' | \\
        sed 's/1|1/1/g' | \\
        sed 's/0|1/NA/g' | \\
        sed 's/1|0/NA/g' | \\
        sed 's/.|./NA/g'  | \\
        sed 's/0\\/0/-1/g' | \\
        sed 's/1\\/1/1/g'  | \\
        sed 's/0\\/1/NA/g' | \\
        sed 's/1\\/0/NA/g' | \\
        sed 's/.\\/./NA/g' > ${group}_${maf}_Genotype_Matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${group}_${maf}_Genotype_Matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}