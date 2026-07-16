---
title: "Bug Report"
date: 2026-03-23
params:
  title: "Incorrect LOCO Marker Count Test Assertions"
  nf_version: ""
  executor: "local"
  command: "TEST_DB_DIR=.../tests/integration_data/test/db TEST_WORK_DIR=.../tests/.nf-work TEST_LEGACY_ASSESSMENT=.../simulation_assessment_results.tsv TEST_DB_ASSESSMENT=.../db_simulation_assessment_results.tsv Rscript tests/run_tests.R"
  error: |
    ── 1. Failure (test-cross_validation.R:349:5): loco GWA has fewer markers than inbred for same marker
    Expected loco n_markers (7016) < inbred n_markers (7016) for population=ce.test.200strains maf=0.05 < `inbred_row$n_markers`.
    Actual comparison: 7016 >= 7016
    Difference: 0 >= 0

    ── 2. Failure (test-db_structure.R:568:5): loco GWA marker count is less than marker set size (GCTA ex
    Expected loco n_markers (7016) < marker set size (7016) for population=ce.test.200strains — GCTA exclusions present < `ms_row$n_markers[1]`.
    Actual comparison: 7016 >= 7016
    Difference: 0 >= 0
  config: ""
  context: |
    The failing tests were written under the assumption that GCTA silently excludes near-singular markers from LOCO GWA output, causing LOCO marker counts to be strictly less than the marker set size (and less than INBRED marker counts). A prior fix corrected the LOCO marker count so that it now matches the full marker set. The tests were never updated to reflect this fix and now incorrectly assert a strict inequality that is no longer valid.

    **Fix:** Remove the two failing assertions (or convert them to equality checks where appropriate) and revise the testing documentation to reflect that LOCO and INBRED marker counts are expected to be equal after the fix.
format: gfm
engine: knitr
---

# Bug Report: Incorrect LOCO Marker Count Test Assertions

**Date:** 2026-03-23
**Nextflow version:** N/A (test-suite only)
**Executor:** local

---

## Description

Two integration tests incorrectly assert that LOCO GWA output should contain *fewer* markers than the inbred mapping and the full marker set. This assumption was valid before a prior fix that corrected the LOCO marker count; the fix made LOCO marker counts equal to the marker set size, which causes both tests to fail.

- `test-cross_validation.R:349` — asserts `loco n_markers < inbred n_markers`
- `test-db_structure.R:568` — asserts `loco n_markers < marker set n_markers`

Both tests now observe equality (`7016 >= 7016`, difference: 0) and fail.

## Command Used

```bash
# Step 1 — collect integration test data
bash tests/collect_test_data.sh

# Step 2 — run integration tests
TEST_DB_DIR=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/integration_data/test/db \
  TEST_WORK_DIR=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/.nf-work \
  TEST_LEGACY_ASSESSMENT=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/integration_data/test/simulation_assessment_results.tsv \
  TEST_DB_ASSESSMENT=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/integration_data/test/db_simulation_assessment_results.tsv \
  Rscript tests/run_tests.R
```

## Error Output

```
── 1. Failure (test-cross_validation.R:349:5): loco GWA has fewer markers than inbred for same marker set
Expected loco n_markers (7016) < inbred n_markers (7016) for population=ce.test.200strains maf=0.05 < `inbred_row$n_markers`.
Actual comparison: 7016 >= 7016
Difference: 0 >= 0

── 2. Failure (test-db_structure.R:568:5): loco GWA marker count is less than marker set size (GCTA exclusions present)
Expected loco n_markers (7016) < marker set size (7016) for population=ce.test.200strains — GCTA exclusions present < `ms_row$n_markers[1]`.
Actual comparison: 7016 >= 7016
Difference: 0 >= 0
```

## Config / Profile

```groovy
// N/A — failure is in the R test suite, not pipeline execution
```

## Additional Context

**Root cause:** The tests encode the old behavior (GCTA silently excludes near-singular markers from LOCO output → LOCO count < marker set count). A fix upstream corrected the LOCO marker count so it now equals the full marker set size. The tests were never updated.

**Fix:** Remove the two failing strict-inequality assertions. If LOCO vs. INBRED count parity is still worth testing, replace with an equality check. Update the testing documentation to state that LOCO and INBRED marker counts are expected to be equal after the fix.

**Affected files:**
- `tests/testthat/test-cross_validation.R` — line 349
- `tests/testthat/test-db_structure.R` — line 568

---

*Generated with the Quarto GitHub Issue Framework. Submit via `gh issue create --title "Fix incorrect LOCO marker count test assertions" --label "bug" --body-file issues/incorrect-loco-marker-count-test/bug_report.md`*
