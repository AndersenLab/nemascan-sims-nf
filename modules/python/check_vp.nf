\
process PYTHON_CHECK_VP {

    label 'python_check_vp'
    tag "${nqtl} ${rep} ${h2} ${effect} ${mode} ${group}_${maf}"

    input:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix)
    tuple path(tmp_pheno_file), path(hsq_file) // From GCTA_MAKE_GRM.out.pheno_and_hsq
    path par_file                              // Original .par file, needs to be piped into this process
    path check_vp_script                       // Path to bin/check_vp.py

    output:
    tuple val(group), val(maf), val(nqtl), val(effect), val(rep), val(h2), val(mode), val(suffix), emit: params
    tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"), path(par_file), emit: pheno
    path "versions.yml", emit: versions

    script:
    def final_pheno_name = "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"
    """
    python3 ${check_vp_script} --check_vp ${hsq_file} --simulated_phenos ${tmp_pheno_file}

    if [ -f new_phenos.temp ]; then
        mv new_phenos.temp ${final_pheno_name}
    else
        cp ${tmp_pheno_file} ${final_pheno_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def final_pheno_name = "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"
    """
    touch ${final_pheno_name}
    // par_file is an input, so it's expected to exist for the stub.
    // No need to touch par_file in stub output as it's passed through.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
} 