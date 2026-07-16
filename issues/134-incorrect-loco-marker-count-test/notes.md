
## Command Run:

Run the following command to collect test data by running the test profile of the pipeline:
```bash
bash tests/collect_test_data.sh
```

Then perform integration testing

```bash
$ TEST_DB_DIR=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/integration_data/test/db \
  TEST_WORK_DIR=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/.nf-work \
  TEST_LEGACY_ASSESSMENT=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/integration_data/test/simulation_assessment_results.tsv \
  TEST_DB_ASSESSMENT=/Users/ryanmckeown/Desktop/fix-merge-unit-with-integration/tests/integration_data/test/db_simulation_assessment_results.tsv \
  Rscript tests/run_tests.R
```

## Error Message:

```txt
Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Some features are not enabled in this build of Arrow. Run `arrow_info()` for more information.
The repository you retrieved Arrow from did not include all of Arrow's features.
You can install a fully-featured version by running:
`install.packages('arrow', repos = 'https://apache.r-universe.dev')`.

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: DBI

Attaching package: ‘readr’

The following objects are masked from ‘package:testthat’:

    edition_get, local_edition


Attaching package: ‘data.table’

The following object is masked from ‘package:purrr’:

    transpose

The following objects are masked from ‘package:dplyr’:

    between, first, last

assessment_cross_validation: ..........................
cross_validation: ....................................1....
cv_pool: SSSSSSSS
db_structure: ..........................................................................................S.2.

══ Skipped ═══════════════════════════════════════════════════════════════════════════════════════════
1. traits/causal_genotypes/ directory exists and contains parquet files (test-cv_pool.R:73:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

2. causal genotype parquets have correct schema (test-cv_pool.R:87:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

3. non-marker causal variants appear as FN rows with NA log10p in assessment (test-cv_pool.R:109:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

4. non-marker FN rows have significant=NA (not FALSE) in assessment (test-cv_pool.R:126:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

5. Simulated.QTL.VarExp is populated for non-marker FN rows (test-cv_pool.R:141:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

6. non-marker causal variants appear in causal_variants/ but not in markers DB view (test-cv_pool.R:157:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

7. designate_qtl() classifies non-marker FN rows as Missed.CV (test-cv_pool.R:188:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

8. calculate_simrep_performance() runs correctly on cv_pool assessment data (test-cv_pool.R:206:3) - Reason: TEST_CV_POOL=true not set — skipping cv_pool integration tests

9. all expected populations appear in mappings_metadata and are queryable (test-db_structure.R:488:5) - Reason: TEST_EXPECTED_POPULATIONS not set — skipping multi-population assertions

══ Failed ════════════════════════════════════════════════════════════════════════════════════════════
── 1. Failure (test-cross_validation.R:349:5): loco GWA has fewer markers than inbred for same marker 
Expected loco n_markers (7016) < inbred n_markers (7016) for population=ce.test.200strains maf=0.05 < `inbred_row$n_markers`.
Actual comparison: 7016 >= 7016
Difference: 0 >= 0

── 2. Failure (test-db_structure.R:568:5): loco GWA marker count is less than marker set size (GCTA ex
Expected loco n_markers (7016) < marker set size (7016) for population=ce.test.200strains — GCTA exclusions present < `ms_row$n_markers[1]`.
Actual comparison: 7016 >= 7016
Difference: 0 >= 0
```

## RCA

The failure here is due to an incorrect test. We made a fix so that local markers have the correct number of markers. So it is incorrect to assume that local GWAS should have fewer markers than inbred for the same marker set. 

## Fix

Remove the failed test, assume and check for local having fewer numbers, and revise the testing documentation to address this.

