# Feature request: pre-transfer manifest + archival workflow

## Summary

Add a reusable tool, `pretransfer_manifest.sh`, that stages a completed
NemaScan-sims-nf output directory on the HPC side as a single reproducible
`tar.zst` archive with a per-file SHA-256 manifest and run-level metadata.
The purpose is to make Rockfish → local-workstation transfers safe (full
integrity verification), efficient (one big file instead of thousands of
small Parquet files), and auditable (reproducible hashes, provenance
metadata) without introducing new runtime dependencies.

## Motivation

A completed production run (see `2026-04-03-production_three_species_sims.qmd`)
produces a ~11 GB output directory containing:

- The Parquet database at `{outputDir}/db/` — five subdirectories of
  small Parquet files, including a Hive-partitioned `mappings/` tree.
  Thousands of files, most well under 1 MB.
- One or two top-level assessment TSVs (`db_simulation_assessment_results.tsv`
  and optionally `simulation_assessment_results.tsv`).

The current transfer procedure, documented at Phase 8 of
`2026-04-03-production_three_species_sims.qmd`, is a naked
`rsync -av --progress` with **no checksum verification** and no record
of what was transferred:

```bash
rsync -av --progress \
  <user>@<host>:/vast/eande106/<user>/nemascan-sims-nf/Sims_.../ \
  Sims_.../
```

Problems with this:

1. **Small-file overhead.** Rsync pays a per-file protocol cost. Thousands
   of small Parquet files are slow to enumerate and transfer one-by-one
   compared to one bulk-read of a single large file.
2. **No integrity guarantee.** Silent corruption during transfer, on disk,
   or after a partial interrupt is invisible. The user cannot prove that
   the local copy is byte-identical to the HPC copy.
3. **No provenance.** Nothing records which pipeline commit produced the
   output, when the transfer was staged, or who ran it.
4. **No atomic point-in-time snapshot.** Files on HPC can change under
   rsync's feet if the source dir is somehow still being written.

## Goals

- **Atomic integrity** — a single hash certifies the entire payload.
- **Bit-level audit** — a per-file manifest lets you point at any single
  corrupted Parquet file after the fact.
- **Reproducibility** — re-running the stage step on unchanged source
  produces the same archive hash. Drift in the archive hash means source
  drift.
- **Resumability** — partial transfers can resume without restarting.
- **Provenance** — every staging run records the pipeline git commit,
  timestamp, host, user, source byte count, and tool versions.
- **No source deletion** — the tool never removes source files or scratch
  data. Cleanup is a separate, manual, post-verification step.

## Non-goals

The following are explicitly **out of scope** for this feature:

- Archival to cloud storage (S3, GCS, Box, Dropbox, etc.)
- Globus integration
- Multi-tier backups / retention policies
- Encryption at rest
- Transfer of work-directory artifacts (intermediate GRM/PLINK files)
- Automation of the transfer itself — the tool only stages and prints
  copy-pasteable instructions; the user runs `rsync` by hand.

These can be layered on later if needed. The minimum viable tool is the
staging script plus a verification protocol that works with stock SSH
and coreutils on both ends.

## User-facing requirements

- **Single invocation.** One bash command on the HPC side takes an output
  directory and produces all four sidecar files in the staging location.
- **Staging size.** Archive is ~7–9 GB for a production run (input ~11 GB
  with zstd -9 on already-Snappy-compressed Parquet). Roughly doubles
  temporary `/vast/` usage during staging.
- **Standard tool chain.** Uses only bash, GNU coreutils, GNU tar (≥1.28),
  zstd, sha256sum — all present on Rockfish login nodes.
- **Copy-pasteable verify steps.** The script prints the exact commands
  the user needs to run on the local side; no ambiguity.
- **Fail loud.** Non-zero exit codes on every error path; `set -euo
  pipefail` throughout; trap-on-error with line numbers.

## Related files

- `pretransfer_manifest.sh` — the tool (repo root).
- `docs/transfer-and-archival.qmd` — user-facing usage page.
- `issues/pretransfer-manifest-archive/design.md` — technical design.
- `issues/pretransfer-manifest-archive/plan.md` — implementation checklist.
- `2026-04-03-production_three_species_sims.qmd` — Phase 8/9 edited to use
  the new workflow.
- `conf/rockfish.config:206-211` — source directory conventions on Rockfish.
