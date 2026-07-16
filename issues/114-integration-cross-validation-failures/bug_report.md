# Bug Report

2026-03-20

# Bug Report: Integration Cross-Validation Failures

**Date:** 2026-03-20 **Nextflow version:** version 24.10.4, build 5934
(20-01-2025 16:47 UTC) **Executor:** local

------------------------------------------------------------------------

## Description

Errors and failures observed during integration tests, specifically in
`test-cross_validation.R` and `test-db-structure.R`.

## Command Used

Run pipeline (requires Docker) and collect outputs

``` bash
bash tests/collect_test_data.sh --clean
```

Then run integration tests with the collected data:

``` bash
TEST_DB_DIR=tests/integration_data/db \
TEST_WORK_DIR=tests/.nf-work \
TEST_LEGACY_ASSESSMENT=tests/integration_data/simulation_assessment_results.tsv \
TEST_DB_ASSESSMENT=tests/integration_data/db_simulation_assessment_results.tsv \
Rscript tests/run_tests.R
```

## Error Output

    assessment_cross_validation: ................
    assessment_var_exp: .................
    assessment: .....................................................................
    cross_validation: 123.....4
    db_structure: .......................5.......................................
    genotype_storage: ...................
    hash_generation: ...........................................................................
    library_check: .................
    lockfile: ...S
    module_sourcing: ...............
    read_bim: ...
    read_raw_gwa: ..........
    safe_log10p: .........
    schema_changes: ......................................
    significance: ...
    trait_storage: ............................................
    write_gwa_to_db: ........................................

    ══ Skipped ══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
    1. renv.lock contains all packages used in the r_packages container scope (test-lockfile.R:30:3) - Reason: renv.lock is a placeholder stub — run scripts/generate_renv_lock.sh to populate it

    ══ Failed ═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
    ── 1. Error (test-cross_validation.R:79:5): DB mapping P and BETA values match raw GWA source files ─────────────────────────────────────────────────────
    Error in `if (!grepl("^[0-9a-f]{20}$", trait_hash)) {     stop("trait_hash must be a 20-character lowercase hex string (pass $hash, not $hash_string)") }`: the condition has length > 1
    Backtrace:
        ▆
     1. └─global generate_mapping_id(params) at test-cross_validation.R:79:5

    ── 2. Error (test-cross_validation.R:139:5): DB SE values match raw GWA source files ────────────────────────────────────────────────────────────────────
    Error in `if (!grepl("^[0-9a-f]{20}$", trait_hash)) {     stop("trait_hash must be a 20-character lowercase hex string (pass $hash, not $hash_string)") }`: the condition has length > 1
    Backtrace:
        ▆
     1. └─global generate_mapping_id(params) at test-cross_validation.R:139:5

    ── 3. Error (test-cross_validation.R:172:5): DB AF1 values match raw GWA source files ───────────────────────────────────────────────────────────────────
    Error in `if (!grepl("^[0-9a-f]{20}$", trait_hash)) {     stop("trait_hash must be a 20-character lowercase hex string (pass $hash, not $hash_string)") }`: the condition has length > 1
    Backtrace:
        ▆
     1. └─global generate_mapping_id(params) at test-cross_validation.R:172:5

    ── 4. Error (test-cross_validation.R:260:3): marker set count matches number of markers in mappings ─────────────────────────────────────────────────────
    Error in `paste0("v=2|population=", population, "|maf=", sprintf("%.10f", as.numeric(maf)), "|species=", species, "|vcf_release_id=", vcf_release_id, "|ms_ld=", sprintf("%.10f", as.numeric(ms_ld)))`: argument "vcf_release_id" is missing, with no default
    Backtrace:
        ▆
     1. └─global read_marker_set(pop, maf_val, db_dir) at test-cross_validation.R:260:3
     2.   └─global get_markers_path(...)
     3.     └─global generate_marker_set_id(...)
     4.       └─base::paste0(...)

    ── 5. Error (test-db_structure.R:143:3): marker set parquet has correct schema ──────────────────────────────────────────────────────────────────────────
    Error in `get_markers_path(population, maf, species, vcf_release_id, ms_ld, base_dir)`: argument "vcf_release_id" is missing, with no default
    Backtrace:
        ▆
     1. └─global read_marker_set(expected_population, 0.05, db_dir) at test-db_structure.R:143:3
     2.   └─global get_markers_path(...)
     3.     └─global generate_marker_set_id(...)
     4.       └─base::paste0(...)

    ══ DONE ═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
    Error in `x$end_reporter()`:
    ! Failures detected.
    Backtrace:
         ▆
      1. └─testthat::test_dir(...)
      2.   └─testthat:::test_files(...)
      3.     └─testthat:::test_files_serial(...)
      4.       └─testthat::with_reporter(...)
      5.         └─reporter$end_reporter()
      6.           └─testthat:::o_apply(self$reporters, "end_reporter")
      7.             └─base::lapply(objects, f)
      8.               └─testthat (local) FUN(X[[i]], ...)
      9.                 └─x$end_reporter(...)
     10.                   └─testthat:::o_apply(self$reporters, "end_reporter")
     11.                     └─base::lapply(objects, f)
     12.                       └─testthat (local) FUN(X[[i]], ...)
     13.                         └─x$end_reporter(...)
     14.                           └─cli::cli_abort("Failures detected.")
     15.                             └─rlang::abort(...)
    Execution halted

## Config / Profile

Integration tests run locally with Docker, using the `test` profile
which includes the `docker` config.

## Additional Context
