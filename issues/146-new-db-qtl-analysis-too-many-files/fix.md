# Fix Plan — P0 structural fix for `DB_MIGRATION_ANALYZE_QTL` FD exhaustion

## Context

See `rca.md` in this directory for the full root-cause analysis. Briefly: every `analyze_qtl` and `assess_sims` R process calls `query_for_threshold_analysis(mapping_id, base_dir)`, which routes through `get_db_connection()` → `open_mapping_db()` and executes:

```r
CREATE VIEW mappings AS
SELECT * FROM read_parquet('.../db/mappings/**/data.parquet',
                           hive_partitioning = true, union_by_name = true)
```

`union_by_name = true` forces DuckDB to open every parquet file in the tree at query-**prepare** time to compute the unified schema. At 16,608 files × ~100 concurrent tasks, the per-task file-descriptor demand exceeds both the per-process `ulimit -n` and the system-wide `fs.file-max` budget on Rockfish compute nodes, and every task in flight crashes with `IO Error: Cannot open file ... Too many open files in system`.

Each `analyze_qtl`/`assess_sims` task already knows its `population` and computes its own `mapping_id` locally — it never needs the global view. This plan replaces the glob-view read with a targeted single-file read of `db/mappings/population={pop}/mapping_id={id}/data.parquet`. Per-task file opens drop from ≈16,608 to ≤3. The database is untouched; the fix is a read-path change only.

## Goal

1. Add a new function `query_mapping_direct()` in `R/queries.R` that returns the same schema as `query_for_threshold_analysis()` but reads a single known partition file directly (no DuckDB view, no glob, no `union_by_name` scan).
2. Swap the call sites in `analyze_qtl.R:108` and `assess_sims.R:203` to use the new function.
3. Leave `query_for_threshold_analysis()` and `open_mapping_db()` intact — they are still used by `R/queries.R:121, 143, 217, 272, 376, 439` (bulk queries, statistics, catalog functions) and by the test suite. Only the two hot-path call sites move.

## Design

### New function: `query_mapping_direct()`

**Location:** `R/queries.R`, inserted after `query_for_threshold_analysis()` (currently ends at line 357).

**Signature:**

```r
query_mapping_direct <- function(mapping_id, population, ms_meta, base_dir = "data/db")
```

- `mapping_id` — 20-char hex string, identifies the single partition to read.
- `population` — population identifier (e.g. `"ce.full"`). Same value the caller already has in `opt$group` / `params$population`.
- `ms_meta` — the named list returned by `read_marker_set_metadata(population, maf, base_dir)`. Both `analyze_qtl.R` (line 78) and `assess_sims.R` (line 83) already load this *before* the hot-path call, so the caller has it in scope with zero extra work. It carries `maf`, `species`, `vcf_release_id`, `ms_ld` — everything needed to locate the markers parquet.
- `base_dir` — database root, same as today.

**Returns:** A dataframe with the exact same columns and column order as `query_for_threshold_analysis()`:

```
marker, mapping_id, P, BETA, SE, var.exp, CHROM, POS, AF1, population, maf
```

ordered by `(CHROM, POS, marker)`. `var.exp` is `NA_real_` (matching the current `NULL AS "var.exp"` behaviour at `R/queries.R:318`).

**Behaviour (step by step):**

1. Construct the partition file path via the existing `get_partition_path(population, mapping_id, base_dir)` helper at `R/database.R:1695-1703`, then append `data.parquet`. Error with a clear message if the file is missing (helps catch typos / half-written DBs early).
2. Read the file with `arrow::read_parquet()`. The raw file holds five columns — `marker, AF1, BETA, SE, P` — per `mappings_schema()` at `R/database.R:266-274`. If the file is empty, return an empty dataframe with the full schema so downstream `nrow()` checks still work.
3. Call `read_marker_set(population, ms_meta$maf, ms_meta$species, ms_meta$vcf_release_id, ms_meta$ms_ld, base_dir)` at `R/database.R:712` to load `marker → (CHROM, POS)` lookup data.
4. `dplyr::left_join` on `marker` to attach `CHROM` and `POS`.
5. `dplyr::mutate` to add `mapping_id`, `population`, `maf = as.numeric(ms_meta$maf)`, and `` `var.exp` = NA_real_ ``.
6. `dplyr::select` to reorder columns exactly as `query_for_threshold_analysis()` returns them.
7. `dplyr::arrange(CHROM, POS, marker)` to match the legacy `ORDER BY` at `R/queries.R:339`.
8. **Preserve the CHROM:POS uniqueness check** currently at `R/queries.R:343-354`. Copy it verbatim into the new function — it guards against LOCO deduplication bugs in the write path and the RCA explicitly called this out as a must-preserve invariant.

**Cost:** at most three file opens — one mapping partition file, one marker set parquet, and one marker set metadata file (already opened by the caller before calling this function). No DuckDB connection, no glob scan, no schema unification.

**Why a new function instead of `direct = TRUE` on `query_for_threshold_analysis()`:**
- `query_for_threshold_analysis()` takes `(mapping_id, base_dir, con)`. The new path needs `population` and `ms_meta` — different input shape.
- Side-by-side functions make the call-site diff explicit and reviewable.
- Both can coexist. `query_for_threshold_analysis()` remains correct (and remains useful for bulk queries that genuinely want a view over many partitions); we just stop calling it from the two hot paths.

### What stays untouched

- **`open_mapping_db()` (`R/queries.R:31-92`)** — still used by `query_db_mappings_summary()`, `query_bulk_for_threshold_analysis()`, and five other places in `R/queries.R`. The glob-view pattern is correct for aggregate queries that really do need to scan multiple partitions; it is only wrong for the per-task `analyze_qtl` / `assess_sims` hot path.
- **`query_for_threshold_analysis()` (`R/queries.R:308-357`)** — kept as-is for any remaining callers (none in the hot path after this change, but the function is referenced by tests and potentially by future bulk-analysis scripts).
- **`get_threshold_params()` (`R/database.R:169-202`)** — does NOT touch `open_mapping_db()`; it already reads marker set metadata directly via `read_marker_set_metadata()`. No change needed.
- **`.db_cache` (`R/database.R:67`)** — unused in the new path, harmless for other consumers. Leave it.
- **Database schema, write path, hash generation, main.nf channels, Nextflow process definitions** — no changes.

## Files to modify

### 1. `R/queries.R`

**Add** a new function `query_mapping_direct()` after line 357 (the closing brace of `query_for_threshold_analysis()`). Approximate size: 60–80 lines including doc comment.

Sketch (not final code):

```r
#' Read a single mapping's data by direct partition file read
#'
#' Returns the same schema as query_for_threshold_analysis(), but reads
#' db/mappings/population={pop}/mapping_id={id}/data.parquet directly via
#' arrow::read_parquet() rather than creating a DuckDB view over the full
#' partition tree. This avoids the `union_by_name = true` schema scan that
#' opens every parquet file in the database on first CREATE VIEW.
#'
#' Use this in per-task hot paths (analyze_qtl.R, assess_sims.R) where the
#' task already knows its population and mapping_id. For aggregate/bulk
#' queries that need multiple mappings, prefer query_for_threshold_analysis()
#' or query_bulk_for_threshold_analysis().
#'
#' @param mapping_id Mapping ID (20-char hex) identifying the target partition
#' @param population Population identifier (e.g. "ce.full")
#' @param ms_meta Named list from read_marker_set_metadata(); must contain
#'   maf, species, vcf_release_id, ms_ld
#' @param base_dir Database root directory
#' @return Dataframe with columns: marker, mapping_id, P, BETA, SE, var.exp,
#'   CHROM, POS, AF1, population, maf (ordered by CHROM, POS, marker)
query_mapping_direct <- function(mapping_id, population, ms_meta, base_dir = "data/db") {
  # 1. Locate the partition file
  partition_path <- get_partition_path(population, mapping_id, base_dir)
  mapping_file <- file.path(partition_path, "data.parquet")
  if (!file.exists(mapping_file)) {
    stop(glue::glue("Mapping partition not found: {mapping_file}"))
  }

  # 2. Read raw mapping rows (5 cols: marker, AF1, BETA, SE, P)
  mapping_rows <- as.data.frame(arrow::read_parquet(mapping_file))

  if (nrow(mapping_rows) == 0) {
    return(data.frame(
      marker = character(), mapping_id = character(),
      P = numeric(), BETA = numeric(), SE = numeric(),
      `var.exp` = numeric(), CHROM = character(), POS = integer(),
      AF1 = numeric(), population = character(), maf = numeric(),
      check.names = FALSE, stringsAsFactors = FALSE
    ))
  }

  # 3. Load marker set for CHROM/POS lookup
  markers <- read_marker_set(
    population,
    as.numeric(ms_meta$maf),
    ms_meta$species,
    ms_meta$vcf_release_id,
    as.numeric(ms_meta$ms_ld),
    base_dir
  )

  # 4. Join + reshape to match query_for_threshold_analysis() output schema
  result <- mapping_rows %>%
    dplyr::left_join(
      markers %>% dplyr::select(marker, CHROM, POS),
      by = "marker"
    ) %>%
    dplyr::mutate(
      mapping_id = mapping_id,
      population = population,
      maf        = as.numeric(ms_meta$maf),
      `var.exp`  = NA_real_
    ) %>%
    dplyr::select(
      marker, mapping_id, P, BETA, SE, `var.exp`,
      CHROM, POS, AF1, population, maf
    ) %>%
    dplyr::arrange(CHROM, POS, marker)

  # 5. Preserve the LOCO dedup guard from query_for_threshold_analysis:343-354
  dup_check <- result %>%
    dplyr::group_by(CHROM, POS) %>%
    dplyr::filter(dplyr::n() > 1)
  if (nrow(dup_check) > 0) {
    n_dups <- nrow(dup_check)
    example <- paste0(dup_check$CHROM[1], ":", dup_check$POS[1])
    stop(glue::glue(
      "CHROM:POS uniqueness violation in mapping {mapping_id}: ",
      "{n_dups} duplicate rows (e.g. {example}). ",
      "Check for LOCO deduplication in write_gwa_to_db.R"
    ))
  }

  result
}
```

Final implementation should match the existing code style in `R/queries.R` (backticks for `var.exp`, lowercase `as.data.frame`, etc.).

### 2. `modules/db_migration/analyze_qtl/resources/usr/bin/analyze_qtl.R`

**Line 108** — replace:

```r
mapping_data <- query_for_threshold_analysis(mapping_id, opt$base_dir)
```

with:

```r
mapping_data <- query_mapping_direct(mapping_id, params$population, ms_meta, opt$base_dir)
```

`params$population` comes from line 66; `ms_meta` from line 78. Both are already in scope — no new reads, no new args, no new library sources.

Nothing else in `analyze_qtl.R` needs to change. `get_threshold_params()` at line 117 does not touch `open_mapping_db()`, so it stays. `nrow(mapping_data)` at line 128 still returns the GCTA output row count (the partition file *is* the GCTA output), preserving the `mckeowr1/issue111` BF-denominator fix.

### 3. `modules/db_migration/assess_sims/resources/usr/bin/assess_sims.R`

**Line 203** — replace:

```r
mapping_data <- query_for_threshold_analysis(mapping_id, opt$base_dir)
```

with:

```r
mapping_data <- query_mapping_direct(mapping_id, params$population, ms_meta, opt$base_dir)
```

Same pattern as `analyze_qtl.R`. `params$population` comes from line 66; `ms_meta` from line 83. The downstream `nrow(mapping_data)` BF-threshold use at line 214 is unchanged.

### 4. No changes to

- `R/database.R`
- `R/analysis.R`, `R/assessment.R`, `R/io.R`, `R/utils.R`
- `main.nf`
- `modules/db_migration/analyze_qtl/main.nf` / `modules/db_migration/assess_sims/main.nf`
- `conf/rockfish.config`
- Container images, lockfiles, Nextflow schema
- Any other `modules/db_migration/*` script

## Verification

### A. Local smoke test (fastest feedback loop)

```bash
nextflow run main.nf -profile test
```

The `test` profile builds a tiny DB (13 strains, chr I/II/V, 2 reps). Expected: both `DB_MIGRATION_ANALYZE_QTL` and `DB_MIGRATION_ASSESS_SIMS` complete; `db_simulation_assessment_results.tsv` is produced. A green run confirms the new code path is wired correctly and returns a schema the downstream code accepts.

### B. Integration test suite (schema + content parity)

```bash
bash tests/collect_test_data.sh           # runs the pipeline + collects fixtures
TEST_DB_DIR=tests/integration_data/db \
TEST_WORK_DIR=tests/.nf-work \
TEST_LEGACY_ASSESSMENT=tests/integration_data/simulation_assessment_results.tsv \
TEST_DB_ASSESSMENT=tests/integration_data/db_simulation_assessment_results.tsv \
Rscript tests/run_tests.R
```

Critical test files that exercise the new path:

- `tests/testthat/test-db_structure.R` — validates the Parquet DB layout after a real pipeline run; will fail loudly if the direct-read helper returns a malformed dataframe.
- `tests/testthat/test-cross_validation.R` — asserts legacy vs DB path parity on mapping data content; end-to-end coverage of `query_mapping_direct()` via `analyze_qtl.R`.
- `tests/testthat/test-assessment_cross_validation.R` — asserts legacy vs DB assessment output parity; covers the `assess_sims.R` call site.

If all three pass with the old path and still pass with the new path, the fix is content-equivalent.

### C. Schema-parity unit test (recommended, small)

Add a focused test in `tests/testthat/test-cross_validation.R` (or a new file) that, given `TEST_DB_DIR`:

1. Picks one `mapping_id` from `mappings_metadata.parquet`.
2. Looks up its `population` and calls `read_marker_set_metadata()` to get `ms_meta`.
3. Calls `query_for_threshold_analysis(mapping_id, TEST_DB_DIR)` → `legacy_df`.
4. Calls `query_mapping_direct(mapping_id, population, ms_meta, TEST_DB_DIR)` → `direct_df`.
5. Asserts `names(legacy_df) == names(direct_df)`, `nrow(legacy_df) == nrow(direct_df)`, and `all.equal(legacy_df, direct_df)` (tolerating column type coercion if needed — e.g., `POS` may come back as integer from one path and double from the other).

This is the fastest way to catch a subtle schema regression before running the full pipeline.

### D. Rockfish resume test

The failing production run's DB is intact. After the fix is merged:

```bash
cd /path/to/failing/run
nextflow run main.nf -profile rockfish -resume \
    --strainfile strains.tsv --vcf 20231213 \
    --output_dir Sims_ce-cb-ct_nqtl1-50_h2grid_50reps
```

Expected: all `DB_MIGRATION_WRITE_*` and `DB_MIGRATION_AGGREGATE_METADATA` tasks are cached; only `DB_MIGRATION_ANALYZE_QTL` and `DB_MIGRATION_ASSESS_SIMS` re-execute. No FD errors. `db_simulation_assessment_results.tsv` is produced.

A good intermediate-confidence signal before the full run: grab one failing work dir, manually `export R_SOURCE_DIR=...`, source the modified R library files, and reproduce the failing `analyze_qtl.R` invocation on the real 16,608-file database. If it returns a populated dataframe without opening more than a handful of files (`lsof -p <pid>` while it runs), the fix is validated at production scale.

## Risks and mitigation

| Risk | Likelihood | Mitigation |
|---|---|---|
| Subtle column type / order drift between old and new paths | Medium | Schema-parity unit test (C) + existing cross-validation tests (B). Run both before merging. |
| Missing partition file that the old glob-view silently skipped | Low | The old view would also crash on a missing file during the `SELECT ... WHERE` (DuckDB would report zero rows). The new path errors earlier with a clearer message, which is an improvement. |
| `ms_meta` schema drift in future (e.g., field rename) breaks the new helper | Low | `ms_meta` is already used identically by `read_genotype_matrix()` and `read_phenotype_data()` in `assess_sims.R:133-143`. Any drift would break those first. |
| LOCO dedup regression | Low | The new helper copies the `CHROM:POS` uniqueness check verbatim from `R/queries.R:343-354`. |
| Some other caller still hits the glob-view under concurrency and triggers the same FD bug | Medium | Other callers (`R/queries.R:121, 143, 217, 272, 376, 439`) are used by bulk/aggregate paths (`query_db_mappings_summary`, `query_bulk_for_threshold_analysis`, `query_population_statistics`, etc.), none of which are called in tight per-task loops from the current `main.nf`. A future P1 stop-gap (remove `union_by_name = true` from `open_mapping_db`) is available if this becomes a problem. Out of scope for this fix. |

## Rollback

Straightforward revert:

1. `git revert` the commit that changed `analyze_qtl.R:108` and `assess_sims.R:203`.
2. Optionally leave `query_mapping_direct()` in `R/queries.R` — it's dead code but harmless, and keeps the helper available for future use.

No database migration, no schema change, no downstream script changes — nothing to undo beyond the two one-line call-site edits.

## Commit / branch plan

**Branch name:** `fix/analyze-qtl-too-many-files` (or similar — matches the issue directory `issues/new-db-qtl-analysis-too-many-files/`).

**Commit sequence** (three small commits for reviewability):

1. `Add query_mapping_direct() for single-partition mapping reads`
   — new function in `R/queries.R` only; no call-site changes yet.
2. `Switch analyze_qtl.R and assess_sims.R to query_mapping_direct()`
   — both call-site changes in a single commit so they land atomically.
3. `Add schema-parity test for query_mapping_direct()` *(optional but recommended)*
   — the unit test from section C.

Include a link to `issues/new-db-qtl-analysis-too-many-files/rca.md` in the PR description so reviewers have the full context without re-deriving it.

## Documentation updates

- **`CLAUDE.md`** — add a bullet under "Development Notes" noting that `analyze_qtl`/`assess_sims` use `query_mapping_direct()` by design to avoid the `union_by_name` FD blow-up; bulk queries may still use `query_for_threshold_analysis()` / `open_mapping_db()`.
- **`docs/database-structure.qmd`** — optional note in the mappings section about the single-partition-read pattern for per-task reads.
- **`docs/changelog.qmd`** — entry describing the fix and referencing the RCA.

No new parameter → `nextflow_schema.json` needs no update.

## Estimated scope

- ~70–90 lines added in `R/queries.R` (new function + doc comment).
- 1 line changed in `analyze_qtl.R`.
- 1 line changed in `assess_sims.R`.
- ~40 lines added in `tests/testthat/` (schema-parity unit test).
- ~10 lines of documentation updates.

Total: small, well-bounded, single-purpose PR.
