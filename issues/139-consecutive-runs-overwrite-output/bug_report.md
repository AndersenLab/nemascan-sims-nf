---
title: "Bug Report"
date: 2026-03-26
params:
  title: "Consecutive same-day pipeline runs silently overwrite prior output directory"
  nf_version: ">=24.10.0"
  executor: "local / rockfish"
  command: "nextflow run main.nf -profile test"
  error: ""
  config: ""
  context: ""
format: gfm
engine: knitr
---

# Bug Report: Consecutive same-day pipeline runs silently overwrite prior output directory

**Date:** 2026-03-26
**Nextflow version:** >=24.10.0
**Executor:** local / rockfish

---

## Description

When the pipeline is run twice on the same calendar day without running
`tests/collect_test_data.sh` (or otherwise moving/archiving the output) between runs, the
second run writes to the same date-stamped output directory
(`Analysis_Results-{YYYY-MM-DD}/`) as the first run. Nextflow publish actions overwrite
existing files silently. The first run's integration test data — including
`db_simulation_assessment_results.tsv` and all Parquet DB files — is lost.

This is especially costly when running with the full CaeNDR VCF (multi-hour runs) and can
cause confusion: the DB on disk reflects the second run, but the researcher may believe they
are examining results from the first run.

## Command Used

```bash
# First run
nextflow run main.nf -profile test

# (No collect_test_data.sh step here — easy to forget)

# Second run same day — silently overwrites Analysis_Results-YYYY-MM-DD/
nextflow run main.nf -profile test
```

## Error Output

```
No error is raised. The pipeline completes successfully and the output directory is
silently overwritten. The first run's results are unrecoverable unless the Nextflow
work directory still contains the original published files.
```

## Config / Profile

```groovy
// nextflow.config — default output directory uses date only (no time component)
// outputDir = "Analysis_Results-${new java.util.Date().format('yyyy-MM-dd')}"
```

## Additional Context

- The default output directory is `Analysis_Results-{YYYY-MM-DD}`. Two runs on the same
  date resolve to the same path.
- The `tests/collect_test_data.sh` workflow moves pipeline outputs into
  `tests/integration_data/` — if this step is skipped, the data in
  `Analysis_Results-{date}/` will be overwritten by the next run.
- The risk is highest on HPC (Rockfish) where runs use the full VCF and take several hours
  to complete. A researcher who queues a second run before archiving the first will lose
  their results.
- `--output_dir` can be passed explicitly to avoid the collision, but there is no guard or
  warning when the target directory already exists.

## Proposed Fix

Add a startup check (in `main.nf` or a validation workflow) that detects when the resolved
`outputDir` already exists and either:

1. **Abort with a clear error** asking the user to pass `--output_dir` with a unique path, or
2. **Append a run index or short timestamp suffix** (e.g., `Analysis_Results-2026-03-26_2`
   or `Analysis_Results-2026-03-26T1430`) to avoid the collision automatically.

A warning in the `collect_test_data.sh` docs and `docs/testing-local-e2e.qmd` noting that
the script must be run before repeating a same-day pipeline run would also reduce confusion.

---

*Submit via `gh issue create --title "Consecutive same-day runs silently overwrite output directory" --label "bug" --body-file issues/consecutive-runs-overwrite-output/bug_report.md`*
