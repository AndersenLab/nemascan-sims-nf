process PLINK_RECODE_VCF {

    label 'plink_recode_vcf'
    tag "${meta.id}_${maf}"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val mito_name
    each val(maf)

    output:
    tuple val(meta.id), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), emit: plink
    tuple val(meta.id), val(maf), path("recoded.vcf.gz"), path("recoded.vcf.gz.tbi")                                                                                                   emit: vcf
    tuple val(meta.id), val(maf), path("markers.txt")                                                                                                                                  emit: markers
    path "versions.yml",                                                                                                                                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cp ${vcf} recoded.vcf.gz
    cp ${vcf_index} recoded.vcf.gz.tbi

    plink --vcf recoded.vcf.gz \\
        --snps-only \\
        --biallelic-only \\
        --maf ${maf} \\
        --set-missing-var-ids @:# \\
        --indep-pairwise 50 10 0.8 \\
        --geno \\
        --not-chr ${mito_name} \\
        --allow-extra-chr

    plink --vcf recoded.vcf.gz \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${maf} \\
        --set-missing-var-ids @:# \\
        --extract plink.prune.in \\
        --geno \\
        --recode \\
        --out TO_SIMS \\
        --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
        sort -k1,1d -k2,2n > markers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$( plink --version |& head -n 1 | cut -f 2' )
    END_VERSIONS
    """

    stub:
    """
    touch TO_SIMS.bed
    touch TO_SIMS.bim
    touch TO_SIMS.fam
    touch TO_SIMS.map
    touch TO_SIMS.nosex
    touch TO_SIMS.ped
    touch TO_SIMS.log
    touch recoded.vcf.gz
    touch recoded.vcf.gz.tbi
    touch markers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$( plink --version |& head -n 1 | cut -f 2' )
    END_VERSIONS
    """
}