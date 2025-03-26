process LOCAL_COMPILE_EIGENS {
    
    label "local_compile_eigens"
    tag "${group} ${maf}"

    input:
    tuple val(group), val(maf), path("*_${group}_${maf}_independent_snvs.csv")

    output:
    tuple val(group), val(maf), path("${group}_${maf}_total_independent_tests.txt"), emit: tests

    script:
    """
    cat *independent_snvs.csv | \\
        grep -v inde | \\
        awk '{s+=\$1}END{print s}' > ${group}_${maf}_total_independent_tests.txt
    """

}
