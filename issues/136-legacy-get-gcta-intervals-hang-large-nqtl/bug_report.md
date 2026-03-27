---
title: "Bug Report"
date: 2026-03-23
---

# Bug Report: `Get_GCTA_Intervals.R` hangs at `readr::write_tsv(processed_mapping)` for large nQTL

**Date:** 2026-03-23
**Nextflow version:** 24.10.4
**Executor:** local (Docker)
**Profile:** `test_variable` (with `--legacy_assess`)
**Related:** Issue #73, Issue #81

---

## Description

`R_GET_GCTA_INTERVALS` hangs indefinitely at `readr::write_tsv(processed_mapping, ...)` for
large-nQTL, low-h² parameter combinations. Observed with `nqtl=30, h2=0.2, inbred, nopca,
EIGEN threshold`. The container runs at ~100% CPU with **0B block I/O** for 20+ minutes —
R is serializing the data frame in memory and never reaches the actual disk write.

The root cause is the size of `processed_mapping`. The `process_mapping_df()` function builds
`correlation_df` by joining the full GWA output (22,582 markers) with the wide-format genotype
matrix (~200 strains per marker), producing a data frame with ~22,582 rows × 200+ columns
(~4.5M cells). `readr::write_tsv` serializes this in-memory before writing, which causes
pathological memory/CPU behavior for large nQTL combinations where all rows are retained.

This is compounded by `dplyr::left_join(..., copy = TRUE)` at line 611 (`Get_GCTA_Intervals.R`),
which may trigger additional data copying behavior.

## Command Used

```bash
bash tests/collect_test_data.sh --profile test_variable
# Which runs:
NXF_VER=24.10.4 nextflow run main.nf -profile test_variable,docker --legacy_assess \
  -work-dir tests/.nf-work
```

## Error Output

No error — the process hangs silently. Last output in `.command.err`:

```
Combining interval positions into a single data frame.

Combined interval positions data frame created with 1 rows.

[1] "Joining Ve correlation data with interval positions..."

[1] "Joined Ve correlation data with interval positions."

[1] "Finished processing mapping data."

[1] "Saving processed mapping data..."
```

The script never prints `"Saved processed mapping data."` (line 663).

Docker container stats while hung:
```
CPU %: 99.70%    MEM: 235.1MiB / 15.6GiB    BLOCK I/O: 0B / 0B
```

## Config / Profile

```groovy
test_variable {
    process.executor = 'local'
    params {
        nqtl   = "data/test/variable_architecture/nqtl.csv"   // includes 30, 50
        h2     = "data/test/variable_architecture/h2.csv"     // includes 0.2
        effect = "data/test/variable_architecture/effect_sizes.csv"
        reps   = 2
    }
}
```

## Additional Context

**Hang location:** `bin/Get_GCTA_Intervals.R` line 658:
```r
readr::write_tsv(processed_mapping,
  file = glue::glue("{trait_name}_{args[12]}_{args[13]}_{args[11]}_processed_{label}_mapping.tsv"),
  col_names = T
)
```

**Why large nQTL is worse:** `process_mapping_df()` constructs `correlation_df` by joining
the mapping data with the genotype matrix (wide format, one column per strain). For nQTL=30
(or nQTL=50) with all 22,582 markers × 200 strains, the resulting object is ~4.5M cells.
The small-nQTL, high-h² combos that ran successfully in the first pipeline run (nQTL=5,
h2=0.8) likely succeeded because R's internal string serialization stays within
practical limits for smaller nQTL cases where fewer markers pass filtering.

**Impact:**
- The `test_variable` profile cannot be run with `--legacy_assess` without hitting this hang
  for nQTL ∈ {15, 20, 30, 50} at low h² values
- The `simulation_assessment_results.tsv` written to `Analysis_Results-*/` is incomplete
  (only covers parameter combos that completed before the hang)
- `test-assessment_cross_validation.R` cannot be fully validated for the `test_variable`
  profile unless the legacy path is repaired or the profile is run with `--no-legacy`

**Workaround:** Run `test_variable` with `--no-legacy`:
```bash
bash tests/collect_test_data.sh --profile test_variable --no-legacy
```
This skips the legacy assessment entirely. The `test-assessment_cross_validation.R` tests
will skip (no `TEST_LEGACY_ASSESSMENT`), but all DB structure, cross-validation, and unit
tests run normally.

---

*Submit via `gh issue create --title "Get_GCTA_Intervals.R hangs at write_tsv for large nQTL (nqtl≥15)" --label "bug" --body-file issues/legacy-get-gcta-intervals-hang-large-nqtl/bug_report.md`*
