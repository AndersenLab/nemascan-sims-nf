process PYTHON_CHECK_VP {

    label 'python_check_vp'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix)
    tuple path("tmp.pheno"), path(hsq_in), path(par_in)
    path check_vp_script

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), emit: params
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path(par_in), emit: pheno
    path "versions.yml", emit: versions

    script:
    def final_pheno_name = "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"
    """
    python3 ${check_vp_script} --check_vp ${hsq_in} --simulated_phenos tmp.pheno

    mv new_phenos.temp ${final_pheno_name}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def final_pheno_name_stub = "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"
    """
    touch ${final_pheno_name_stub}
    touch ${par_in} // par_in is the path object from input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
