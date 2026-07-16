# Bug Report

2026-03-20

# Bug Report: RF integration missing eigen outputs

**Date:** 2026-03-20 **Nextflow version:** 24.10.4 **Executor:**
Rockfish SLURM (pipeline run) / local macOS (test run)

------------------------------------------------------------------------

## Description

The EIGEN threshold cross-validation test
(`test-cross_validation.R:229`) fails with “EIGEN file not found” when
run locally against an `.nf-work` directory synced from Rockfish via
`rsync`.

**Root cause:** On Rockfish, Nextflow stages the EIGEN file
(`*_total_independent_tests.txt`) into every downstream process work
directory (e.g., `DB_MIGRATION_WRITE_MARKER_SET`, `GCTA_PERFORM_GWA`,
`R_FIND_GENOTYPE_MATRIX_EIGEN`) as an **absolute symlink** pointing back
to the canonical `LOCAL_COMPILE_EIGENS` work directory on `/scratch4/`.
There is only one real copy of the file; all other appearances are
symlinks.

When `rsync -av` copies the `.nf-work` tree locally, it preserves
symlinks as-is. The symlink targets (e.g.,
`/scratch4/eande106/Ryan/nf-work-test/d4/3d34704b.../ce.test.200strains_0.05_total_independent_tests.txt`)
do not exist on the local machine, creating **broken symlinks**.

`list.files(recursive = TRUE)` in `test-cross_validation.R` returns
these broken symlink paths. R’s `file.exists()` returns `FALSE` for
broken symlinks, so `read_eigen_file()` throws an error before the test
can iterate to the one real copy of the file.

**Observed state in local `.nf-work`:** - 27 broken symlinks →
`/scratch4/eande106/.../d4/3d34704b.../ce.test.200strains_0.05_total_independent_tests.txt` -
1 real file at
`tests/.nf-work/d4/3d34704b.../ce.test.200strains_0.05_total_independent_tests.txt`

The test iteration hits a broken symlink before reaching the real file
and errors rather than continuing.

## Command Used

``` bash
# Pipeline run (Rockfish)
NXF_VER=24.10.4 nextflow run main.nf \
  -profile test,rockfish \
  --legacy_assess \
  -work-dir /scratch4/eande106/Ryan/nf-work-test

# Sync (local Mac)
rsync -av --delete \
  login.rockfish.jhu.edu:/scratch4/eande106/Ryan/nf-work-test/ \
  tests/.nf-work/

# Test run (local Mac)
TEST_DB_DIR=tests/integration_data/db \
TEST_WORK_DIR=tests/.nf-work \
TEST_LEGACY_ASSESSMENT=tests/integration_data/simulation_assessment_results.tsv \
TEST_DB_ASSESSMENT=tests/integration_data/db_simulation_assessment_results.tsv \
Rscript tests/run_tests.R 2>&1 | tee tests/integration_results_20260320.log
```

## Error Output

    ── 1. Error ('test-cross_validation.R:229:5'): EIGEN threshold from DB matches pipeline EIGEN value
    Error in `read_eigen_file(eigen_file)`:
      EIGEN file not found: /Users/ryanmckeown/Desktop/nemascan-sims-nf/tests/.nf-work/0c/e3dc0c8a2a1c4fc17377fc8201070a/ce.test.200strains_0.05_total_independent_tests.txt
    Backtrace:
        ▆
     1. └─global read_eigen_file(eigen_file) at test-cross_validation.R:229:5

The path exists as a broken symlink:

    lrwxrwxrwx tests/.nf-work/0c/e3dc0c8a.../ce.test.200strains_0.05_total_independent_tests.txt
      -> /scratch4/eande106/Ryan/nf-work-test/d4/3d34704b.../ce.test.200strains_0.05_total_independent_tests.txt

## Config / Profile

``` groovy
// conf/rockfish.config — relevant stageIn behavior (Nextflow default)
// Nextflow stages input files as symlinks by default (stageInMode = 'symlink')
// On shared POSIX filesystems this creates absolute symlinks to the canonical
// work directory, which break when the work tree is relocated to another host.
```

## Additional Context

- The pipeline itself ran successfully on Rockfish; all SLURM jobs
  completed.
- The EIGEN value IS present in `marker_set_metadata.parquet` (DB path
  is unaffected).
- The one real file
  (`d4/3d34704b.../ce.test.200strains_0.05_total_independent_tests.txt`)
  was correctly rsync’d; the test would pass if the iteration reached it
  first.
- This is an HPC-to-local transfer problem, not a pipeline correctness
  issue.

**Candidate fixes (two layers):**

1.  **Test robustness (primary):** In `test-cross_validation.R`, skip
    broken symlinks inside the `for (eigen_file in eigen_files)` loop by
    adding `if (!file.exists(eigen_file)) next` before calling
    `read_eigen_file()`. This lets the test find and use the real copy.

2.  **Sync robustness (secondary):** Add `--copy-unsafe-links` to the
    rsync command so that symlinks whose targets lie outside the
    transferred tree (i.e., absolute HPC paths) are replaced with real
    file copies locally. This makes the synced `.nf-work`
    self-contained.

Both fixes should be applied: the test fix is defensive against any
future relocation scenario; the rsync fix prevents a confusing state
where the work tree appears populated but files are unreadable.

------------------------------------------------------------------------

*Generated with the Quarto GitHub Issue Framework. Submit via
`gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
