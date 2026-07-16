# Root Cause Analysis — `DB_MIGRATION_ANALYZE_QTL` "Too many open files"

## Summary

On the production Rockfish run `Sims_ce-cb-ct_nqtl1-50_h2grid_50reps`, every upstream database write completed successfully, but 2,024 `DB_MIGRATION_ANALYZE_QTL` tasks failed at query-**prepare** time with `IO Error: Cannot open file ... Too many open files in system`. The root cause is structural: `open_mapping_db()` registers a DuckDB view over **every** partition in `db/mappings/**/data.parquet` with `union_by_name = true`, forcing DuckDB to open all 16,608 parquet files on the first query prepare inside each R process — regardless of which single `mapping_id` the task actually needs. Replicated across ~100 concurrent R processes from the `array = 100` SLURM job array, the file-descriptor demand collides with the per-process / per-node limit and every task in flight crashes. The fix is to stop asking DuckDB to scan the whole tree: each `analyze_qtl` / `assess_sims` task already knows its `population` and `mapping_id` at launch time and can read the single partition file directly.

## Failure evidence

### From `issues/new-db-qtl-analysis-too-many-files/bug_report.md`

- **Run:** `Sims_ce-cb-ct_nqtl1-50_h2grid_50reps` on Rockfish
- **Failing process:** `DB_MIGRATION_ANALYZE_QTL`
- **Failure count:** 2,024 tasks
- **Upstream state:** all writes completed; the Parquet DB is intact

### From a failing work directory

Work dir: `/scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct3/59/75afe787ceda437d78566953c36b6c`

```
2026-04-09 01:26:37 [INFO] Analyzing mapping: 2e59b3684db2b4a2e256 with threshold: BF
Error: rapi_prepare: Failed to prepare query   CREATE VIEW mappings AS
  SELECT * FROM read_parquet('/vast/.../db/mappings/**/data.parquet',
                             hive_partitioning = true, union_by_name = true)
Error: IO Error: Cannot open file ".../mappings/population=ce.full/mapping_id=4012d3f36a69c83d744a/data.parquet":
  Too many open files in system
```

**This task was analyzing mapping `2e59b3684db2b4a2e256`, but crashed trying to open `4012d3f36a69c83d744a` — a completely unrelated partition.** That cross-contamination is the smoking gun: DuckDB is opening files the task does not need.

### Database scale (runtime verification)

```bash
$ find .../db/mappings -name 'data.parquet' | wc -l
16608

$ find .../db/mappings -mindepth 1 -maxdepth 1 -type d
.../db/mappings/population=cb.full
.../db/mappings/population=ce.full
.../db/mappings/population=ct.full
```

**16,608 parquet files across three populations** (`cb.full`, `ce.full`, `ct.full`). Every `analyze_qtl` task attempts to open all of them.

## Mechanism

Six-step chain from the task entry point to the file-descriptor blow-up. All line numbers are against the current tree.

1. **Per-task entry.** Each `DB_MIGRATION_ANALYZE_QTL` process launches a fresh `Rscript analyze_qtl.R` and calls `query_for_threshold_analysis(mapping_id, base_dir)` at `modules/db_migration/analyze_qtl/resources/usr/bin/analyze_qtl.R:108`. The script only needs rows for one known `mapping_id` (computed locally at line 104).

2. **Cached DuckDB connection.** That call routes through `get_db_connection(base_dir, use_cache = TRUE)` at `R/queries.R:311`, which on a fresh R process calls `open_mapping_db()` at `R/database.R:78-103`.

3. **Global view over the whole partition tree.** `open_mapping_db()` at `R/queries.R:55-67` executes:

   ```r
   glob_pattern <- file.path(mappings_dir, "**", "data.parquet")
   DBI::dbExecute(con, glue::glue("
     CREATE VIEW mappings AS
     SELECT * FROM read_parquet('{glob_pattern}', hive_partitioning = true, union_by_name = true)
   "))
   ```

   There is no `WHERE` clause and no per-task path scoping — the view is over the entire mappings tree.

4. **`union_by_name = true` forces an eager schema scan.** DuckDB must read the parquet footer of **every** matching file to compute the unified column list before the `CREATE VIEW` statement can be prepared. The R error message (`rapi_prepare: Failed to prepare query CREATE VIEW`) confirms the crash happens during `DBI::dbPrepare` / `dbExecute`, not during the later `SELECT ... WHERE mapping_id = ...` query. Each opened footer consumes one OS file descriptor until the prepare returns.

5. **The cached connection holds descriptors.** `R/database.R:78-103` caches the connection in `.db_cache`. DuckDB keeps its file references alive through the duration of the R process, so the task that opens the descriptors is the task that holds them.

6. **Parallelism multiplies the demand past the kernel limit.** `conf/rockfish.config:167-174` gives the `db_migration_analyze_qtl` label `cpus = 1`, `array = 100`, and **no** `maxForks`. The global `executor.queueSize = 100` at `conf/rockfish.config:200-203` caps Nextflow submission but not concurrent execution. With ~100 R processes each opening 16,608 files, the aggregate system descriptor demand is on the order of **1.66 million open FDs**, which collides with the per-process `ulimit -n` (typically 1,024 / 4,096) on the first few dozen opens and with the system-wide `fs.file-max` / `fs.nr_open` shortly after. DuckDB reports both EMFILE and ENFILE as "Too many open files in system".

**The `WHERE m.mapping_id = '{mapping_id}'` clause at `R/queries.R:338` never enters the picture.** It runs only after the view is prepared, so hive partition pruning cannot undo the schema-unification file opens at prepare time.

## Scale math

| Quantity | Value |
|---|---|
| Parquet files in `db/mappings/` | **16,608** (verified) |
| Populations | 3 (`cb.full`, `ce.full`, `ct.full`) |
| Concurrent `analyze_qtl` tasks per SLURM array | ≤100 (`array = 100`) |
| File opens per task (at CREATE VIEW prepare) | ≈16,608 |
| Peak system-wide open FDs | ≈1.66 M |
| Typical Linux `ulimit -n` (soft) | 1,024 – 4,096 |
| Typical Linux `/proc/sys/fs/file-max` | 0.8 M – 2 M |

Per-process exhaustion is hit almost immediately (16,608 ≫ 1,024). Even if `ulimit -n` were raised to 65,536, the *system-wide* FD budget would still be saturated by a handful of concurrent tasks on the same compute node.

The exhaustion is deterministic, not probabilistic — any run whose `db/mappings/` tree grows past ~1,000 files will start to fail under concurrent read, and the failure rate scales with both database size and job concurrency.

## Why upstream processes succeeded

Write-path modules (`DB_MIGRATION_WRITE_MARKER_SET`, `DB_MIGRATION_WRITE_GWA_TO_DB`, `DB_MIGRATION_WRITE_TRAIT_DATA`, `DB_MIGRATION_AGGREGATE_METADATA`) each open only their own output file plus a small fixed set of inputs. Per-task FD cost is **O(1)**, not **O(N_mappings)**, so they scale cleanly no matter how large the database has already grown. The issue is specific to the first read module that constructs the global DuckDB view — and that's `analyze_qtl`.

## Why `DB_MIGRATION_ASSESS_SIMS` is the next shoe to drop

`modules/db_migration/assess_sims/resources/usr/bin/assess_sims.R:203` calls `query_for_threshold_analysis()` through the same cached-connection / global-view path. If `analyze_qtl` were simply retried or skipped, `assess_sims` would hit the identical failure mode on the same input scale. **Any fix must land in both scripts (or in the shared R library function they both call).**

## Hypothesis falsification — what this is *not*

- **Not a corrupted parquet file.** The error is `Cannot open file: Too many open files`, not a parse error. The failing file opens fine in isolation (confirmed by the fact that upstream `WRITE_GWA_TO_DB` wrote it successfully).
- **Not a database-barrier timing bug.** `DB_MIGRATION_AGGREGATE_METADATA.out.summary` at `main.nf:702` correctly gates `analyze_qtl`; the FD error is downstream of a fully-populated DB.
- **Not a DuckDB version regression.** The view pattern predates the failing run. What changed is the database *scale*, not the DuckDB call.
- **Not the `normalizePath("~")` warning.** That's a cosmetic R warning from the container missing a `HOME` dir; it's benign and unrelated to the fatal error.
- **Not per-task randomness.** The failing task was analyzing `2e59b3684db2b4a2e256` but crashed opening `4012d3f36a69c83d744a` — a completely different mapping, confirming the view is globbing everything rather than scoping to the task's own partition.

## Fix directions (prioritized)

### P0 — Structural fix: single-file direct read by known partition path **(recommended)**

Bypass the DuckDB glob-view entirely for `analyze_qtl` and `assess_sims`. Each task already knows its `population` (`opt$group`) and computes its own `mapping_id` (`analyze_qtl.R:104`), which uniquely identifies a file at `db/mappings/population={pop}/mapping_id={id}/data.parquet`. The read becomes three small, known-path operations:

1. `arrow::read_parquet()` on the single partition file — 5 columns (`marker`, `AF1`, `BETA`, `SE`, `P`).
2. Reuse existing `read_marker_set_metadata()` + `read_marker_set()` in `R/database.R` to attach `CHROM`, `POS`, `maf`, and the `marker_set_id` (needed by the `nrow(mapping_data)` BF-denominator fix from `mckeowr1/issue111`).
3. Re-attach `mapping_id` and `population` locally — they're known at call time.

**Files to change (deferred to a follow-up implementation plan):**

- `R/database.R` — add `read_mapping_by_id_direct(population, mapping_id, base_dir)` near `get_partition_path()` at approximately line 1703. Reuse `get_partition_path()`, `read_marker_set_metadata()`, `read_marker_set()`.
- `R/queries.R` — add an opt-in `direct = TRUE` path on `query_for_threshold_analysis()` (or a sibling function) that returns the same schema but reads via the new helper. Preserve the `CHROM:POS` uniqueness check at `R/queries.R:343-354` and the `NULL AS "var.exp"` default at `R/queries.R:318`.
- `modules/db_migration/analyze_qtl/resources/usr/bin/analyze_qtl.R:108` — call the direct path.
- `modules/db_migration/assess_sims/resources/usr/bin/assess_sims.R:203` — same change.

**Effect.** Per-task file opens drop from ≈16,608 to ≤5 (one mapping file + one marker file + one marker-set-metadata file + `mappings_metadata.parquet`). The `union_by_name` schema scan disappears from the hot path. Scales to any future database size without re-tuning.

**Gotchas to preserve during implementation:**

- The raw partition file does **not** contain `mapping_id`, `population`, `CHROM`, `POS` — all must be re-attached (as they already are via the current LEFT JOINs in `query_for_threshold_analysis`).
- The `nrow(mapping_data)` BF-denominator fix from `mckeowr1/issue111` must continue to work (the direct read returns the same row count as the GCTA output, since the partition file *is* the GCTA output).
- `safe_log10p()` in `R/analysis.R` is still computed downstream from `P` — unchanged.
- `open_mapping_db()` is still used by other `R/queries.R` functions that genuinely need a cross-mapping aggregate view (e.g., `query_bulk_for_threshold_analysis`, `mapping_statistics`). Don't remove it wholesale — just opt the hot-path tasks out.

### P1 — Structural fix lite: drop `union_by_name = true`

Removing `union_by_name = true` from `R/queries.R:64` (and from the `marker_set_metadata` view at `R/queries.R:85`) allows DuckDB to infer the schema from a single file and defer the rest to query-time scans with hive partition pruning. This is a smaller diff than P0 and would meaningfully reduce the file-open count, but it is a **stop-gap**:

- Any future schema drift across partitions silently changes behavior.
- The scan still touches more than the one file the task needs.
- It does not remove the "share one view across unrelated queries" footgun.

Consider P1 only if P0 is blocked for some reason.

### P2 — Tactical band-aid: cap concurrency

Add `maxForks = 10` (or similar) to the `db_migration_analyze_qtl` and `db_migration_assess_sims` labels at `conf/rockfish.config:167-174` and `:175-182`. This does not fix the per-task FD blow-up — a single task with ≈16,608 open files still exceeds a 1,024 per-process limit — but it reduces the *system-wide* pressure so that small-to-medium databases (where per-process stays under `ulimit -n`) don't trip the node limit. **Label explicitly as a band-aid.** It does not help the current 16,608-file failure and should not be relied on past a few thousand partitions.

### P3 — Environment tuning: raise `ulimit -n` in the container / SLURM job

Raising `ulimit -n` via `singularity.runOptions` or `SLURM --propagate=NOFILE` might paper over small runs, but it's fragile:

- Depends on kernel `fs.nr_open` cap.
- Adds 10×–100× more descriptors per task regardless of whether they're actually needed.
- Does nothing for system-wide `fs.file-max` exhaustion.

Only useful as defence-in-depth *after* P0 is shipped, not as the primary fix.

### Not recommended

- **DuckDB file-handle cache tuning.** DuckDB does not expose a knob to limit how many files the `union_by_name` schema scan opens concurrently. The only reliable way to reduce it is to not ask for the scan.
- **Rerunning the pipeline with `-resume` and no code change.** The same code path will fail again at the same scale. The DB is already intact — do not re-run upstream work.

## Recommended next step

Open a follow-up implementation plan for **P0**:

1. Add `read_mapping_by_id_direct()` to `R/database.R` (or a new helper colocated with `query_for_threshold_analysis`).
2. Switch `analyze_qtl.R:108` and `assess_sims.R:203` to the direct path.
3. Local `-profile test` smoke test to confirm schema parity with the current path.
4. Small cross-validation check: run both the old (glob-view) and new (direct) path on a local test DB and `all.equal()` the returned frames for a handful of mapping IDs.
5. Resume the failed Rockfish run with `-resume` once the fix is in — all upstream work is preserved in the existing Parquet DB, so only `DB_MIGRATION_ANALYZE_QTL` and `DB_MIGRATION_ASSESS_SIMS` need to re-execute.

No database regeneration is required; the fix is a read-path change only.
