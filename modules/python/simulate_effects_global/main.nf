process PYTHON_SIMULATE_EFFECTS_GLOBAL {

  label 'python_simulate_effects_global'
  tag "${nqtl} ${rep} ${effect}${group}_${maf}"

  input:
  tuple val(group), val(maf),
        path("TO_SIMS.bed"),    path("TO_SIMS.bim"),     path("TO_SIMS.fam"),
        path("TO_SIMS.map"),    path("TO_SIMS.nosex"),   path("TO_SIMS.ped"),   path("TO_SIMS.log"),
        path(gm), path(n_indep_tests),
        path("CV_TO_SIMS.bed"), path("CV_TO_SIMS.bim"),  path("CV_TO_SIMS.fam"),
        path("CV_TO_SIMS.map"), path("CV_TO_SIMS.nosex"), path("CV_TO_SIMS.ped"), path("CV_TO_SIMS.log")
  path create_causal_qtls
  each rep
  each nqtl
  each effect

  output:
  tuple val(group), val(maf), val(nqtl), val(effect), val(rep),
        path("causal.variants.sim.${nqtl}.${rep}.txt"),
        emit: causal
  tuple val(group), val(maf),
        path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"),
        path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests),
        emit: plink
  tuple val(group), val(maf),
        path("CV_TO_SIMS.bed"),  path("CV_TO_SIMS.bim"),  path("CV_TO_SIMS.fam"),
        path("CV_TO_SIMS.map"),  path("CV_TO_SIMS.nosex"), path("CV_TO_SIMS.ped"), path("CV_TO_SIMS.log"),
        emit: cv_plink
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  """
      python ${create_causal_qtls} CV_TO_SIMS.bim ${nqtl} ${effect}
      mv causal_vars.txt causal.variants.sim.${nqtl}.${rep}.txt

  cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version 2>&1 | cut -d' ' -f2 )
          numpy: \$( python -c "import numpy; print(numpy.__version__)" )
          pandas: \$( python -c "import pandas; print(pandas.__version__)" )
  END_VERSIONS
  """

  stub:
  """
    touch causal.variants.sim.${nqtl}.${rep}.txt
    touch CV_TO_SIMS.bed
    touch CV_TO_SIMS.bim
    touch CV_TO_SIMS.fam
    touch CV_TO_SIMS.map
    touch CV_TO_SIMS.nosex
    touch CV_TO_SIMS.ped
    touch CV_TO_SIMS.log

  cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: stub
          numpy: stub
          pandas: stub
  END_VERSIONS

  """
}
