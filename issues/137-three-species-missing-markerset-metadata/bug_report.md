# Bug Report

2026-03-27

# Bug Report: Concurrent write race in marker_set_metadata.parquet loses population records on HPC

**Date:** 2026-03-27 **Nextflow version:** 24.10.x **Executor:** slurm
(Rockfish)

------------------------------------------------------------------------

## Description

During HPC testing with the `test_hpc_three_species` profile on
Rockfish, two DB migration processes failed for a subset of populations:

| Process | Failed Populations | Status |
|----|----|----|
| `DB_MIGRATION_WRITE_GWA_TO_DB` | `ct.hpc.popB` (inbred pca, inbred nopca) | FAILED (exit 1), cascade ABORTED 20+ tasks |
| `DB_MIGRATION_WRITE_TRAIT_DATA` | `ct.hpc.popB`, `ce.hpc.popB` | FAILED (exit 1), retried 3× each |

All C. briggsae populations and the `.hpc` / `.hpc.popA` groups for C.
elegans and C. tropicalis completed successfully. Only the `.popB`
subpopulations for CE and CT failed.

The error in both processes is identical:

    Marker set metadata not found for population='<pop>', maf=0.05

## Command Used

``` bash
nextflow run main.nf -profile test_hpc_three_species,rockfish -work-dir /scratch4/eande106/Ryan/nf-work-test-hpc-three-species
```

## Error Output

### `DB_MIGRATION_WRITE_GWA_TO_DB` (ct.hpc.popB, inbred pca)

    Error: Marker set metadata not found for population='ct.hpc.popB', maf=0.05
    in /vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Analysis_Results-20260326/db.
    Ensure DB_MIGRATION_WRITE_MARKER_SET completed before resuming.
    If VCF/species parameters changed, a full re-run (not -resume) is required.

### `DB_MIGRATION_WRITE_TRAIT_DATA` (ct.hpc.popB)

Same error message, same root cause. The R script `write_trait_data.R`
calls `read_marker_set_metadata("ct.hpc.popB", 0.05)` which returns
`NULL` because the population’s record is absent from
`marker_set_metadata.parquet`.

## Config / Profile

``` groovy
// rockfish.config — DB migration process labels
withLabel: db_migration_write_marker_set {
    time = "30.minute"
    cpus = 2
    memory = "8G"
    errorStrategy = 'retry'
    maxRetries = 3
    array = 100
}
```

### Three-species strainfile structure

The `test_hpc_three_species` profile uses
`data/test/hpc_strains_three_species.txt` with 9 population groups (3
species × {full, popA, popB}):

| Group                            | Species      |
|----------------------------------|--------------|
| ce.hpc, ce.hpc.popA, ce.hpc.popB | c_elegans    |
| cb.hpc, cb.hpc.popA, cb.hpc.popB | c_briggsae   |
| ct.hpc, ct.hpc.popA, ct.hpc.popB | c_tropicalis |

All groups share the same MAF threshold (0.05), producing **9 concurrent
`DB_MIGRATION_WRITE_MARKER_SET` SLURM jobs**.

## Additional Context

Three-species HPC test profile (9 population groups × 1 MAF), SLURM
executor, fresh run (not -resume)

------------------------------------------------------------------------

# Root Cause Analysis

## Summary

**Concurrent read-modify-write race condition** in
`write_marker_set_metadata()` (`R/database.R:1201-1244`). Multiple
`DB_MIGRATION_WRITE_MARKER_SET` SLURM jobs write to the same
`marker_set_metadata.parquet` file simultaneously, causing **lost
updates** — some population records are silently overwritten and never
appear in the final metadata file.

## Mechanism

`write_marker_set_metadata()` uses a non-atomic read-modify-write
pattern:

``` r
# R/database.R lines 1231-1240
if (file.exists(metadata_path)) {
    existing <- as.data.frame(arrow::read_parquet(metadata_path))
    existing_filtered <- existing[existing$marker_set_id != marker_set_id, , drop = FALSE]
    result <- dplyr::bind_rows(existing_filtered, new_record)
} else {
    result <- new_record
}
arrow::write_parquet(result, metadata_path)
```

All 9 `DB_MIGRATION_WRITE_MARKER_SET` tasks on the SLURM cluster target
the same shared file at `{db_output}/marker_set_metadata.parquet`. The
race window is between the `read_parquet()` and `write_parquet()` calls.

### Timeline (illustrative)

    Time    Task A (ce.hpc)                Task B (ce.hpc.popA)           Task C (ce.hpc.popB)
    ────    ───────────────                ────────────────────           ────────────────────
    T1      read metadata → (empty)
    T2                                     read metadata → (empty)
    T3      write [ce.hpc]
    T4                                     write [ce.hpc.popA]            read → [ce.hpc.popA]
                                           ▲ OVERWRITES ce.hpc!
    T5                                                                    write [ce.hpc.popA, ce.hpc.popB]
                                                                          ▲ ce.hpc lost

With 9 concurrent jobs, the final file contains only the records that
survived the last batch of overlapping writes. The “losers” —
populations whose records were read by another task before being
written, then overwritten — are silently dropped.

## Why the barrier doesn’t prevent this

The `ch_marker_barrier` in `main.nf` (lines 590-598) correctly waits for
all `WRITE_MARKER_SET.out.done` signals. However, the barrier only
proves that each task **exited successfully** (emit `val true`), not
that its metadata record **persists in the file**. A task can:

1.  Write its record to the Parquet file
2.  Emit `done` and exit 0
3.  Have its record overwritten by a concurrent task that read the file
    *before* step 1

The barrier guarantees temporal ordering (all marker set tasks finish
before downstream tasks start) but cannot guarantee data integrity when
the write target is a shared mutable file.

## Why this only manifests on HPC

| Environment | Populations | Concurrent marker set tasks | Race probability |
|----|----|----|----|
| `test` (local) | 1 | 1 | None |
| `test_variable` (local) | 1 | 1 | None |
| `test_hpc_three_species` (SLURM) | **9** | **Up to 9** | **High** |

On the local executor, Nextflow schedules tasks sequentially or with
limited parallelism (bounded by `cpus`), and all tasks share a single
process. On SLURM, each task is an independent job writing to the same
file on the `/vast` parallel filesystem (GPFS) with no coordination.

## Evidence

### Failure pattern matches race condition

- All 9 `WRITE_MARKER_SET` tasks completed (the barrier resolved —
  downstream tasks ran), but 2 populations are missing from the metadata
  file
- The missing populations (`ce.hpc.popB`, `ct.hpc.popB`) are the
  “losers” whose records were overwritten by later-finishing concurrent
  writes
- C. briggsae populations all survived — likely because their tasks
  happened to write last (or had less overlap)

### Downstream cascade

- `ct.hpc.popB` failures in `WRITE_GWA_TO_DB` caused Nextflow’s
  `errorStrategy='retry'` to re-attempt (and fail again — metadata still
  missing), then ABORT all dependent downstream tasks
- `ce.hpc.popB` and `ct.hpc.popB` failures in `WRITE_TRAIT_DATA` show
  the same pattern: 3 retries, all fail with the same error

## Questions for investigation

1.  **Confirm marker set task status**: Can you run
    `nextflow log crazy_lalande -f name,status,exit,workdir | grep -i "WRITE_MARKER_SET"`
    to verify all 9 `DB_MIGRATION_WRITE_MARKER_SET` tasks show
    COMPLETED/exit 0?

2.  **Inspect the metadata file**: Can you read the surviving
    `marker_set_metadata.parquet` to see which populations are present?

    ``` r
    arrow::read_parquet(
      "/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Analysis_Results-20260326/db/marker_set_metadata.parquet"
    ) |> dplyr::select(population, maf, marker_set_id)
    ```

    If fewer than 9 rows appear, or if `ce.hpc.popB` / `ct.hpc.popB` are
    absent, this confirms the race condition.

3.  **Check marker data files**: The per-population marker data files
    (`{pop}_{maf}_markers.parquet`) are written to separate files and
    should all exist — can you verify all 9 are present in
    `marker_sets/`?

    ``` bash
    ls -la Analysis_Results-20260326/db/marker_sets/*_markers.parquet | wc -l
    ```

------------------------------------------------------------------------

# Proposed Fix

## Approach: Per-population metadata files

Eliminate the shared mutable file by writing one metadata file per
population/MAF pair, mirroring the existing pattern for marker data
files.

### Current (broken)

All tasks append to a single shared file:

    marker_sets/
      marker_set_metadata.parquet          ← 9 tasks race to read-modify-write
      ce.hpc_0.05_markers.parquet          ← 1 task, no contention
      ce.hpc.popA_0.05_markers.parquet     ← 1 task, no contention
      ...

### Proposed (safe)

Each task writes its own isolated metadata file:

    marker_sets/
      ce.hpc_0.05_metadata.parquet         ← 1 task, no contention
      ce.hpc_0.05_markers.parquet
      ce.hpc.popA_0.05_metadata.parquet    ← 1 task, no contention
      ce.hpc.popA_0.05_markers.parquet
      ...

### Changes required

| File | Change |
|----|----|
| `R/database.R` | `write_marker_set_metadata()` — write to `{pop}_{maf}_metadata.parquet` (no read-modify-write) |
| `R/database.R` | `read_marker_set_metadata()` — read from `{pop}_{maf}_metadata.parquet` directly |
| `R/database.R` | `get_all_marker_set_metadata()` — glob `*_metadata.parquet` and `bind_rows()` |
| `R/database.R` | `get_marker_set_metadata_path()` — accept `population` + `maf` params, or add new per-file path helper |
| `R/queries.R` | Update any DuckDB view that reads `marker_set_metadata.parquet` by name |
| Tests | Update `test-db_structure.R` and `test-schema_changes.R` for new file layout |

### Why this approach

- **Zero race window**: Each task writes only its own file — no shared
  resource
- **No locking needed**: File locks are unreliable on parallel
  filesystems (GPFS/Lustre)
- **No serialization overhead**: All tasks can still run in parallel
- **Consistent with existing patterns**: Marker data
  (`{pop}_{maf}_markers.parquet`) and genotype data
  (`{pop}_{maf}_genotypes.parquet`) already use per-population files
- **Backward compatible**: `get_all_marker_set_metadata()` can fall back
  to reading the single-file format if `*_metadata.parquet` files don’t
  exist

### Alternative considered: `maxForks 1`

Adding `maxForks 1` to the `DB_MIGRATION_WRITE_MARKER_SET` process would
serialize execution and eliminate the race. This is simpler but: -
Reduces parallelism (minor impact with 9 short tasks) - Doesn’t fix the
underlying design flaw — any future concurrent access (e.g., `-resume`
with partial completion) could still trigger the race - Couples
correctness to scheduling, which is fragile

------------------------------------------------------------------------

*Generated with the Quarto GitHub Issue Framework. Submit via
`gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
