---
title: "Bug Report"
date: 2026-03-23
---

# Bug Report: `assess_sims.R` crashes when no QTLs are detected (empty `qtl_regions` file)

**Date:** 2026-03-23
**Nextflow version:** 24.10.4 (via `NXF_VER=24.10.4` in `collect_test_data.sh`)
**Executor:** local
**Profile:** `test_variable`

---

## Description

`DB_MIGRATION_ASSESS_SIMS` crashes when a mapping produces no significant QTLs. This is a
valid and expected outcome (e.g., `nqtl=1, h2=0.4, gamma` effect — low power), but the
pipeline treats it as a fatal error and aborts all downstream assessment.

The failure chain:

1. `analyze_qtl.R` calls `extract_qtl_regions(processed)`.
2. `extract_qtl_regions()` returns `data.frame()` (zero columns, zero rows) when `peak_id`
   is all-NA (no significant markers).
3. `write.table(data.frame(), file, col.names = TRUE)` writes a single newline byte (`0x0a`)
   — column headers cannot be written when the data frame has no columns.
4. `assess_sims.R` calls `data.table::fread(opt$qtl_regions, header = TRUE)` on this
   1-byte file and throws:

   ```
   Error in data.table::fread(opt$qtl_regions, header = TRUE) :
     Input is either empty, fully whitespace, or skip has been set after the last non-whitespace.
   ```

Note: `compile_full_assessment()` in `R/assessment.R` (line 377) already handles
`nrow(qtl_regions) == 0` gracefully — it constructs an empty `peak_info` data frame and
continues. The bug is entirely in the file I/O handoff between `analyze_qtl.R` and
`assess_sims.R`.

**Affected process:** `DB_MIGRATION_ASSESS_SIMS (BF 1 2 0.4 gamma inbred pca ce.test.200strains_0.05)`

## Command Used

```bash
bash tests/collect_test_data.sh --profile test_variable
# Which runs:
NXF_VER=24.10.4 nextflow run main.nf -profile test_variable,docker --legacy_assess \
  -work-dir tests/.nf-work
```

## Error Output

From `.command.err` in work dir `d0/6b529be0c8b4fc5db0e046a837a3b3`:

```
2026-03-23 21:24:41 [INFO] Assessing mapping: 4889f3980364287d4d1e
Error in data.table::fread(opt$qtl_regions, header = TRUE) :
  Input is either empty, fully whitespace, or skip has been set after the last non-whitespace.
Calls: %>% -> as.data.frame -> <Anonymous>
Execution halted
```

Confirmed the `qtl_regions` file is 1 byte (just `0x0a`):
```
$ xxd 1_2_0.4_0.05_gamma_ce.test.200strains_inbred_pca_BF_qtl_regions.tsv
00000000: 0a
```

## Config / Profile

```groovy
test_variable {
    process.executor = 'local'
    params {
        strainfile = "${projectDir}/data/test/test_strains.txt"
        nqtl       = "${projectDir}/data/test/variable_architecture/nqtl.csv"
        h2         = "${projectDir}/data/test/variable_architecture/h2.csv"
        effect     = "${projectDir}/data/test/variable_architecture/effect_sizes.csv"
        reps       = 2
    }
}
```

## Additional Context

**Root cause:** `write.table(data.frame(), file, col.names = TRUE)` produces a 1-byte file
(just a newline) because a zero-column data frame has no column names to emit.
`data.table::fread` requires at least a header row.

**Fix location:** `assess_sims.R` line 114. Add a file-size guard before calling `fread`:

```r
# Before:
qtl_regions <- data.table::fread(opt$qtl_regions, header = TRUE) %>% as.data.frame()

# After:
qtl_regions <- if (file.size(opt$qtl_regions) <= 1L) {
  data.frame()
} else {
  data.table::fread(opt$qtl_regions, header = TRUE) %>% as.data.frame()
}
```

`compile_full_assessment()` already handles `nrow(qtl_regions) == 0` (constructs empty
`peak_info`, returns all-Missed.CV assessment rows). No other changes needed.

The `test_variable` profile is the first profile to exercise low-power parameter
combinations (e.g., `nqtl=1, h2=0.2`) where zero detections are expected. The default
`test` profile uses `nqtl=5, h2=0.8` (high power) and would rarely/never trigger this.

---

*Submit via `gh issue create --title "assess_sims.R crashes on empty qtl_regions file (no QTLs detected)" --label "bug" --body-file issues/assess-sims-empty-qtl-regions-crash/bug_report.md`*
