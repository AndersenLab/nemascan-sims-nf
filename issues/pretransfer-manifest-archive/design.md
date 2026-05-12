# Design: pre-transfer manifest + archival workflow

## Architecture

Three phases, two machines, four sidecar files.

```
[ HPC / Rockfish ]                [ Transit ]           [ Local ]
pretransfer_manifest.sh RUN   →   rsync -P   →   verify + extract + verify
  ├─ RUN.sha256 (per-file)                          ├─ sha256sum -c archive
  ├─ RUN.tar.zst                                    ├─ tar -I zstd -xvf
  ├─ RUN.tar.zst.sha256                             └─ sha256sum -c RUN.sha256
  └─ RUN.manifest.txt (metadata)
```

**Phase 1 — stage on HPC.** A single invocation of the script walks the
output directory, writes a per-file SHA-256 manifest, builds a reproducible
`tar.zst` archive, hashes the archive, and emits a human-readable metadata
file. Source is never touched.

**Phase 2 — transfer.** The user runs `rsync -avP --partial` from the local
machine pulling the four sidecar files. `--partial` preserves interrupted
byte streams so dropped SSH sessions resume from where they left off.

**Phase 3 — verify and extract.** The user runs three commands locally:
archive hash check → extract → per-file hash check. Only after every line
prints `OK` is the local copy trustable. Only then is scratch cleanup on
HPC safe.

## File inventory

For an output directory named `<RUN>` (e.g. `Sims_ce-cb-ct_nqtl1-50_h2grid_50reps`),
the script emits four files into the staging directory:

| File | Contents | Purpose |
|------|----------|---------|
| `<RUN>.tar.zst` | Compressed archive of `<RUN>/` | The payload. |
| `<RUN>.tar.zst.sha256` | Single-line `sha256sum` output for the archive | Fast atomic integrity check — reject bad transfers before extraction. |
| `<RUN>.sha256` | `sha256sum` of every file inside `<RUN>/`, paths relative to `<RUN>`'s parent directory | Per-file audit. After extraction, `sha256sum -c` validates every file. |
| `<RUN>.manifest.txt` | `key: value` metadata (see below) | Provenance and diagnostics. |

The metadata file records:

```
output_dir:       <path>
created_at:       <UTC ISO8601>
hostname:         <hostname>
user:             <user>
pipeline_commit:  <git rev-parse HEAD of the repo containing the script>
file_count:       <N>
source_bytes:     <N>
archive_path:     <path>
archive_bytes:    <N>
archive_sha256:   <hex>
manifest:         <path>
tar_version:      <first line of tar --version>
zstd_version:     <first line of zstd --version>
```

## Why tar + zstd

**Bundling.** The fundamental throughput issue with the pipeline output is
the sheer number of small Parquet files (thousands, most well under 1 MB).
Rsync and SCP both pay a protocol cost per file; on a fast WAN this is
often the dominant term. A single large file reads straight from disk into
the network socket and back out the other side.

**Compression level.** Parquet files use Snappy compression, so further
compression of the Parquet payload is marginal — the DuckDB-emitted Parquet
is already close to the knee of the size curve. Real compression gains
come from the TSV assessment files (plain text, compresses ~5×) and from
tar headers (thousands of ~512-byte blocks per file compress down to
almost nothing). `zstd -9` is the sweet spot:

- Level 3 (default) — ~2 GB savings vs. raw tar, ~5 s/GB
- Level 9 (chosen) — ~2.5 GB savings vs. raw tar, ~15 s/GB, multi-threaded
- Level 19 — ~2.7 GB savings vs. raw tar, ~3 min/GB, diminishing returns

`-T0` uses all available cores. Staging a ~11 GB run takes a few minutes
on a Rockfish login node.

**Atomic integrity.** A single SHA-256 over the archive means corruption
detection is one cheap command. Detection is far easier than recovery,
and we lean into that.

## Reproducibility

The archive flags are chosen so the same source byte-for-byte produces
the same archive bytes:

```
tar --sort=name
    --owner=0 --group=0 --numeric-owner
    --mtime='UTC 2020-01-01'
    -cf -
```

- `--sort=name` — GNU tar ≥ 1.28 feature. Without this, file order depends
  on how `find` walks the directory, which depends on filesystem inode
  order (not stable across runs, not stable across filesystems).
- `--owner=0 --group=0 --numeric-owner` — strips the UID/GID of whoever
  ran the script. Without this, two users staging the same source produce
  different archive bytes.
- `--mtime='UTC 2020-01-01'` — strips per-file mtimes. Without this, a
  file re-touched by the pipeline (e.g. `-resume`) would change the
  archive hash even if contents are unchanged.

Reproducibility is a test oracle: if two stages of unchanged source
produce different archive hashes, something drifted under us.

## Integrity model

Two layers:

1. **Archive SHA-256** (fast-reject) — catches every transit-level
   problem: truncated transfers, single-bit flips, swapped files, wrong
   archive moved into place.
2. **Per-file SHA-256 manifest** (precise-audit) — catches every storage
   problem: bit rot on the local disk, partial extraction, accidental
   edits of extracted files. Critically, it tells you *which* file is
   bad, so you can re-transfer just that one.

Either layer alone is insufficient. The archive hash tells you the transit
was clean but not that the extraction or long-term storage is. The
per-file manifest tells you the on-disk state is good but not that the
bytes you have are the bytes the HPC had (it is generated *from* the local
bytes, not transported independently with the archive).

Both together give the full guarantee. That is why both sidecars ship
together.

## Safety rails

- **No deletes.** The script never removes source files, scratch data, or
  previously-staged artifacts. Every cleanup step is manual and documented.
- **Idempotency guard.** If all four staging artifacts exist and are newer
  than the newest file in source, the script prints "staging is current"
  and exits 0 without work. `--force` bypasses the guard.
- **Fail loud.** `set -euo pipefail`, explicit exit codes (1/2/3/4), `trap`
  prints the failing line. Partial artifacts are left in place for
  inspection on error.
- **Validation first.** Before any work, the script checks that the output
  directory exists, contains a `db/` subdir (sanity check for a pipeline
  output), and that all required tools (`tar`, `zstd`, `sha256sum`, etc.)
  are on `PATH`.
- **GNU tar hard requirement.** Reproducibility depends on GNU tar flags
  that BSD tar does not support. The script detects GNU tar explicitly
  and fails with install instructions rather than silently producing
  non-reproducible archives.

## Provenance cross-check

The Parquet database already records a `strainfile_hash` column in
`marker_set_metadata.parquet` (see `R/database.R` for the
`marker_set_metadata_schema()` function). This is a SHA-256 of the full
strainfile at pipeline submission time, stored inside the DB. The script
does **not** duplicate this; it trusts the DB's internal provenance.

What the script adds is *run-level* provenance at the archive level:

- `pipeline_commit` — the git HEAD of the repo containing the script at
  staging time. Lets you answer "which code produced this archive" later.
- `created_at` — UTC timestamp of the staging run.
- `hostname` + `user` — who staged it and where.
- `tar_version` / `zstd_version` — the tools used, for future decompressors.

Combined, these give a full audit trail from "bytes on local disk" back to
"the exact git commit and the exact strainfile that produced them."

## Alternatives considered

### A. Rsync the directory tree directly with a sidecar per-file manifest

Pros: files remain browsable mid-transfer; resumable mid-file; simpler
mental model.

Cons: **slow**. Rsync's per-file protocol cost dominates for thousands
of small files. No atomic integrity guarantee (the manifest and the files
are independent, so a race during transfer can produce a manifest that
doesn't match what was actually copied). No bit-reproducibility.

Rejected because the production run's file-count dominates the transfer
time budget.

### B. Streaming `tar | ssh "tar -x"` pipeline

Pros: no local disk used on HPC side (no staging); transfer happens in
one shot.

Cons: not resumable (a dropped SSH session means restarting from byte 0);
no checksum verification possible until after extraction; no way to
separate the "build archive" step from the "transfer" step (so a slow
network interacts badly with a slow compressor). No atomic integrity
artifact — nothing to re-verify later.

Rejected because resumability and post-transfer verification are design
requirements.

### C. Globus / Aspera / commercial transfer tools

Pros: purpose-built for large-scale HPC data movement; parallel streams;
built-in integrity.

Cons: requires endpoint setup on both sides (not already available here);
external dependency; out of scope for a simple "snapshot and verify"
tool; users would still want a manifest for later audits.

Can be layered on later if data volume grows beyond what rsync-over-SSH
can handle, but not needed for an 11 GB payload.

### Chosen: staged tar.zst + rsync

Best tradeoff for the user's chosen destination (local workstation) and
packaging (single archive): fast transfer (single large file), resumable
(`rsync --partial`), atomic integrity (archive hash), bit-reproducible
(tar flags), no new dependencies (stock coreutils + GNU tar + zstd).
