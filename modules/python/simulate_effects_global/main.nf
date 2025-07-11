process PYTHON_SIMULATE_EFFECTS_GLOBAL {

  label 'python_simulate_effects_global'
  tag "${nqtl} ${rep} ${effect}${group}_${maf}"

  input:
  tuple val(group), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests)
  path create_causal_qtls
  each rep
  each nqtl
  each effect

  output:
  tuple val(group), val(maf), val(nqtl), val(effect), val(rep), path("causal.variants.sim.${nqtl}.${rep}.txt"), emit: causal
  tuple val(group), val(maf), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log"), path(gm), path(n_indep_tests), emit: plink
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  """
      python ${create_causal_qtls} TO_SIMS.bim ${nqtl} ${effect}
      mv causal_vars.txt causal.variants.sim.${nqtl}.${rep}.txt
    
  cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version 2>&1 | cut -d' ' -f2 )
  END_VERSIONS
  """

  stub:
  """
    touch causal.variants.sim.${nqtl}.${rep}.txt

  cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version 2>&1 | cut -d' ' -f2 )
  END_VERSIONS

  """
}
