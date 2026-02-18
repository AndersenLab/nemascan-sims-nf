process R_FIND_GENOTYPE_MATRIX_EIGEN {

    label 'r_find_genotype_matrix_eigen'
    tag "${chrom} ${group}_${maf}"

    input:
    tuple val(group), val(maf), path(genomatrix)
    path get_genomatrix_eigen
    each chrom

    output:
    tuple val(group), val(maf), val(chrom), path("${chrom}_${group}_${maf}_independent_snvs.csv"), emit: eigen
    path "versions.yml",                                                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat ${genomatrix} | \\
        awk -v chrom="${chrom}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${chrom}_gm.tsv

    Rscript --vanilla ${get_genomatrix_eigen} ${chrom}_gm.tsv ${chrom}

    mv ${chrom}_independent_snvs.csv ${chrom}_${group}_${maf}_independent_snvs.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    touch ${chrom}_gm.tsv
    touch ${chrom}_${group}_${maf}_independent_snvs.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( Rscript --version |& cut -f 4 )
    END_VERSIONS
    """
}