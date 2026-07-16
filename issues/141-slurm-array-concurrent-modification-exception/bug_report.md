---
title: "Bug Report"
date: 2026-04-03
params:
  title: "ConcurrentModificationException in SLURM job array task submission with Singularity"
  nf_version: "24.10.4"
  executor: "slurm (Rockfish)"
  command: "nextflow run main.nf -profile rockfish --strainfile data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/strains_three_species.tsv --vcf 20231213"
  error: ""
  config: ""
  context: ""
format: gfm
engine: knitr
---

# Bug Report: ConcurrentModificationException in SLURM job array task submission with Singularity

**Date:** 2026-04-03
**Nextflow version:** 24.10.4
**Executor:** slurm (Rockfish)

---

## Description

A large production run (3 species, 50 nQTL values, 50 replicates, 4 h2 values) crashes with
`java.util.ConcurrentModificationException` when Nextflow attempts to submit a
`PLINK_UPDATE_BY_H2` task. The error occurs at the **task submission** level inside
Nextflow's SLURM job array collector and Singularity container builder — not during
task execution. Because it is a submission-level exception, `errorStrategy = 'retry'`
does not catch it and the entire session aborts.

The root cause is a thread-safety race between two concurrent consumers of
`GCTA_SIMULATE_PHENOTYPES` output channels. Nextflow DSL2 implicitly broadcasts process
outputs to multiple consumers, but the underlying `ArrayList` backing each tuple may be
shared by reference. When the `TaskArrayCollector` iterates the input file list to build
Singularity `--bind` mount arguments for one consumer, while another consumer's operator
thread concurrently accesses the same list, the iteration fails.

The race is latent — it existed before the triggering commit — and only surfaces under
high concurrency (thousands of GCTA tasks completing in rapid bursts on SLURM).

## Command Used

```bash
#| eval: false
nextflow run main.nf \
  -profile rockfish \
  --strainfile data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/strains_three_species.tsv \
  --nqtl data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/nqtl.csv \
  --h2 data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/h2.csv \
  --effect data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/effect_sizes.csv \
  --reps 50 \
  --cv_maf 0.05 \
  --cv_ld 0.8 \
  --output_dir Sims_ce-cb-ct_nqtl1-50_h2grid_50reps \
  -work-dir /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct
```

## Error Output

```
Apr-03 21:01:56.237 [Actor Thread 26] ERROR nextflow.processor.TaskProcessor -
  Error executing process > 'PLINK_UPDATE_BY_H2 (5 36 0.2 gamma ce.full_0.05)'

Caused by:
  java.util.ConcurrentModificationException

java.util.ConcurrentModificationException: null
    at java.base/java.util.ArrayList$Itr.checkForComodification(ArrayList.java:1013)
    at java.base/java.util.ArrayList$Itr.next(ArrayList.java:967)
    at nextflow.executor.BashWrapperBuilder.createContainerBuilder(BashWrapperBuilder.groovy:706)
    at nextflow.executor.BashWrapperBuilder.makeBinding(BashWrapperBuilder.groovy:286)
    at nextflow.executor.BashWrapperBuilder.buildNew0(BashWrapperBuilder.groovy:410)
    at nextflow.executor.BashWrapperBuilder.buildNew0(BashWrapperBuilder.groovy:172)
    at nextflow.executor.BashWrapperBuilder.build(BashWrapperBuilder.groovy:425)
    at nextflow.executor.GridTaskHandler.prepareLauncher(GridTaskHandler.groovy:107)
    at nextflow.processor.TaskArrayCollector.createTaskArray(TaskArrayCollector.groovy:141)
    at nextflow.processor.TaskArrayCollector.collect(TaskArrayCollector.groovy:103)
    at nextflow.processor.TaskProcessor.submitTask(TaskProcessor.groovy:2339)
    at nextflow.processor.TaskProcessor.checkCachedOrLaunchTask(TaskProcessor.groovy:839)
    at groovyx.gpars.dataflow.operator.ForkingDataflowOperatorActor$1.run(
         ForkingDataflowOperatorActor.java:58)
    at java.base/java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1136)
    at java.base/java.lang.Thread.run(Thread.java:840)

Apr-03 21:01:56.440 [Actor Thread 26] DEBUG nextflow.Session -
  Session aborted -- Cause: java.util.ConcurrentModificationException
```

## Config / Profile

```groovy
// conf/rockfish.config — relevant sections
process {
    executor = "slurm"
    clusterOptions = "-A eande106 -e errlog.txt -N 1"

    withLabel: plink_update_by_h2 {
        time = "10.minute"
        cpus = 4
        memory = "20G"
        errorStrategy = 'retry'
        maxRetries = 3
        array = 100    // SLURM job array batching — triggers the race
    }
}

executor {
    queueSize = 100
    submitRateLimit = 10
}

singularity {
    enabled = true
    autoMounts = true
    runOptions = "--bind ${projectDir}"
}
```

## Additional Context

### Dual consumption pattern (root cause)

`GCTA_SIMULATE_PHENOTYPES` emits three output channels (`.out.params`, `.out.plink`,
`.out.pheno`). Each was consumed by two operators without explicit forking:

1. A `.merge()` chain building `ch_trait` for DB trait data writes (lines 452-467)
2. `PLINK_UPDATE_BY_H2` process invocation (lines 474-478)

Nextflow DSL2 implicitly broadcasts process outputs, but may share the underlying
mutable `ArrayList` by reference. When `TaskArrayCollector` collects PLINK_UPDATE_BY_H2
tasks into a SLURM job array and calls `BashWrapperBuilder.createContainerBuilder` to
iterate input file paths for Singularity `--bind` mounts, the merge chain's operator
thread can concurrently access the same `ArrayList`, causing the
`ConcurrentModificationException`.

The same pattern existed for `GCTA_PERFORM_GWA` outputs (consumed by DB write path, DB
analysis path, and optionally the legacy intervals path).

### Why `errorStrategy = 'retry'` does not help

The exception occurs inside `TaskProcessor.submitTask` → `TaskArrayCollector.collect` —
before the task is submitted to SLURM. Nextflow's retry mechanism only handles
execution-level failures (non-zero exit codes from completed jobs). A submission-level
exception aborts the entire session.

### Triggering conditions

- **High concurrency:** Thousands of GCTA_SIMULATE_PHENOTYPES tasks completing in rapid
  bursts (from a `-resume` of the workflow), flooding the PLINK_UPDATE_BY_H2 submission queue
- **SLURM job arrays:** `array = 100` enables batching, which means `TaskArrayCollector`
  iterates multiple tasks' input lists in sequence while other threads continue processing
- **Singularity:** `BashWrapperBuilder.createContainerBuilder` iterates input file paths
  to build `--bind` mount arguments — this is the specific iteration that races

The race is timing-dependent and may not reproduce on every run, but becomes more likely with larger runs or with runs where many GCTA tasks complete simultaneously like in a `-resume` scenario.

## Fix Applied

Added explicit `.tap{}` channel forks before dual consumption sites, giving each consumer
its own copy of the tuple data:

```groovy
// GCTA_SIMULATE_PHENOTYPES forks
GCTA_SIMULATE_PHENOTYPES.out.params
    .tap { ch_gcta_params_for_trait }
    .set { ch_gcta_params_for_plink }

GCTA_SIMULATE_PHENOTYPES.out.pheno
    .tap { ch_gcta_pheno_for_trait }
    .set { ch_gcta_pheno_for_plink }

GCTA_SIMULATE_PHENOTYPES.out.plink
    .tap { ch_gcta_plink_for_trait }
    .set { ch_gcta_plink_for_plink }
```

The same pattern was applied to `GCTA_PERFORM_GWA` outputs (`.out.params`, `.out.gwa`,
`.out.pheno`), which had 2-3 consumers depending on `--legacy_assess`.

---

*Submit via `gh issue create --title "ConcurrentModificationException in SLURM job array task submission with Singularity" --label "bug" --body-file issues/140-slurm-array-concurrent-modification-exception/bug_report.md`*
