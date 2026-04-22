---
title: "Bug Report"
date: today
params:
  title: "Untitled Bug"
  nf_version: ""
  executor: "local"
  command: ""
  error: ""
  config: ""
  context: ""
format: gfm
engine: knitr
---

# Bug Report: Rename Chromosome Process Termination

---

## Description

The BCF Tools Rename Chromosome Process failed for one replicate. It says the process might have been terminated by the external system.

## Command Used

```bash
nextflow run main.nf \
  -profile rockfish \
  --strainfile data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/strains_panelsize_ce.tsv \
  --nqtl data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/nqtl.csv \
  --h2 data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/h2.csv \
  --effect data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/effect_sizes.csv \
  --reps 50 \
  --cv_maf 0.05 \
  --cv_ld 0.8 \
  --output_dir Sims_panelsize_ce_nqtl5_h2-0.8_50reps \
  -work-dir /scratch4/eande106/Ryan_tmp/nf-work-panelsize
```

## Error Output


### Caused by

```
Caused by:
  Process `BCFTOOLS_RENAME_CHROMS (ce.n100.r05)` terminated for an unknown reason -- Likely it has been terminated by the external system it has been terminated by the external system
```

## Additional Context


### Failed Task Information

* Task ID: `258`
* tag: `ce.n100.r05`

The only task metrics that are reported in the _report.html file are:

* duration: 9m 24s
* realtime: 8m 29s

### Task Resource Allocation 

There is no resource allocation for the `BCFTOOLS_RENAME_CHROMS` process in the configuration file `conf/rockfish.config` 

And no resources are specified within the [process definition]('modules/bcftools/rename_chroms/main.nf')

#### Checking the `.nextlfow/` work dir metadata

From the failed workdir

```bash
cd /scratch4/eande106/Ryan_tmp/nf-work-panelsize/78/4d80229a4907e29340d49b7ecdb34c
```

```bash
$ cat .command.run | grep "#SBATCH"
```

```bash
#SBATCH -J nf-BCFTOOLS_RENAME_CHROMS_(ce.n100.r05)
#SBATCH -o /scratch4/eande106/Ryan_tmp/nf-work-panelsize/78/4d80229a4907e29340d49b7ecdb34c/.command.log
#SBATCH --no-requeue
#SBATCH --signal B:USR2@30
#SBATCH -t 01:00:00
#SBATCH --mem 4096M
#SBATCH -A eande106 -e errlog.txt -N 1
```

I beleive that this resouce allocation matches the default in `conf/rockfish.config`


```conf/rockfish.config
process {
    executor = "slurm"
    clusterOptions = "-A eande106 -e errlog.txt -N 1"
    time = "1.hour"
    cpus = 1
    memory = "4G"
    partition = "parallel"
...
}
```



---

## Root Cause Analysis

### Approach

The Nextflow report HTML (`20260422-10-54-26_report.html`) was parsed to extract task-level metrics for all 149 `BCFTOOLS_RENAME_CHROMS` tasks in the run. The following Python was used to pull status, duration, peak RSS, and error action for each task:

```python
import re

with open('20260422-10-54-26_report.html') as f:
    content = f.read()

pattern = (
    r'\"task_id\":\"(\d+)\"[^{]*\"process\":\"BCFTOOLS_RENAME_CHROMS\"'
    r'[^{]*\"tag\":\"([^\"]+)\"[^{]*\"status\":\"([^\"]+)\"'
    r'[^{]*\"duration\":\"([^\"]+)\"[^{]*\"realtime\":\"([^\"]+)\"'
    r'[^{]*\"peak_rss\":\"([^\"]+)\"[^{]*\"memory\":\"([^\"]+)\"'
    r'[^{]*\"error_action\":\"([^\"]+)\"'
)
matches = re.findall(pattern, content)
```

The failed task (task 258) was also inspected directly to recover its full JSON record.

### Key Results

**Resource usage across all successful tasks:**

| Metric | Successful tasks (n=98) | Allocated |
|--------|------------------------|-----------|
| Peak RSS | ~22–23 MB | 4096 MB |
| Runtime | 140–320 s | 3600 s |

**Failed task (task 258, `ce.n100.r05`) record:**

```json
{
  "task_id": "258",
  "status": "FAILED",
  "exit": "-",
  "duration": "564735",
  "realtime": "509996",
  "peak_rss": "-",
  "peak_vmem": "-",
  "cpus": "1",
  "memory": "4294967296",
  "time": "3600000",
  "error_action": "TERMINATE"
}
```

All resource metrics are `"-"` (not captured), and there is no exit code. This is the signature of a job killed by SLURM at the infrastructure level (node failure or preemption) rather than a process that ran out of memory or time.

**Cascade effect:** Because `BCFTOOLS_RENAME_CHROMS` had no `withLabel` entry in `conf/rockfish.config`, it inherited the default `errorStrategy`, which is `terminate`. When task 258 was killed, Nextflow immediately stopped the pipeline — causing **50+ downstream tasks to be `ABORTED`**.

### Root Cause

The failure was caused by a **transient SLURM infrastructure event** (node failure or preemption), not insufficient resource allocation. The allocated 4 GB and 1-hour limit are both appropriate given that successful tasks used only ~23 MB peak RSS and completed in 2–5 minutes.

The root cause of the pipeline-level failure is the **missing `errorStrategy = 'retry'`** for the `bcftools_rename_chroms` label. Every other labeled process in `conf/rockfish.config` has `errorStrategy = 'retry'` / `maxRetries = 3`, but `BCFTOOLS_RENAME_CHROMS` had no label override, so a single transient kill event terminated the entire run.

### Fix Applied

Added a `withLabel: bcftools_rename_chroms` block to `conf/rockfish.config`:

```groovy
withLabel: bcftools_rename_chroms {
    errorStrategy = 'retry'
    maxRetries = 3
}
```

No resource changes were made — the default 4 GB / 1 CPU / 1 hour are sufficient.

---

*Generated with the Quarto GitHub Issue Framework. Submit via `gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
