#!/usr/bin/env bash
# pretransfer_manifest.sh — stage a NemaScan-sims-nf pipeline output directory
# as a reproducible tar.zst plus a per-file SHA-256 manifest, an archive hash,
# and a metadata sidecar. Operator guide: 2026-04-03-production_three_species_sims.qmd.

set -euo pipefail

[[ $# -lt 1 ]] && { echo "usage: pretransfer_manifest.sh <out_dir> [staging_dir]" >&2; exit 1; }

out_dir="${1%/}"
staging_dir="${2:-./staging}"
staging_dir="${staging_dir%/}"

if [[ ! -d "$out_dir/db" ]]; then
    echo "ERR: ${out_dir} is missing a db/ subdir — not a pipeline output" >&2
    exit 1
fi
mkdir -p "$staging_dir"

RUN="$(basename "$out_dir")"
parent="$(cd "$(dirname "$out_dir")" && pwd)"
RUN_TAR="${staging_dir}/${RUN}.tar.zst"
RUN_TAR_SHA="${staging_dir}/${RUN}.tar.zst.sha256"
RUN_SHA="${staging_dir}/${RUN}.sha256"
RUN_META="${staging_dir}/${RUN}.manifest.txt"

# Disk-space precheck: crude upper bound (source bytes vs free bytes in staging).
# Archive is smaller post-zstd, but this catches "staging filled up mid-zstd".
source_bytes="$(du -sb "$out_dir" | awk '{print $1}')"
free_bytes="$(df --output=avail -B1 "$staging_dir" | tail -n 1 | tr -d ' ')"
if (( free_bytes < source_bytes )); then
    echo "ERR: staging dir ${staging_dir} has ${free_bytes} bytes free, need >= ${source_bytes}" >&2
    exit 1
fi

# Per-file manifest. cd to parent so paths are <RUN>/... — matching how tar
# names entries below. LC_ALL=C is the precondition for stable byte-order sort
# across locales (and for archive reproducibility).
echo "Generating per-file manifest..."
( cd "$parent" && find "$RUN" -type f | LC_ALL=C sort | xargs sha256sum ) > "$RUN_SHA"

# Reproducible tar.zst. The tar flags are the test-oracle reference for
# byte-identical archives across runs; do not change without updating the
# reproducibility check in the runbook.
echo "Creating reproducible tar.zst archive..."
tar --sort=name --owner=0 --group=0 --numeric-owner \
    --mtime='UTC 2020-01-01' -C "$parent" -cf - "$RUN" \
  | zstd -T0 -9 -o "$RUN_TAR"

echo "Hashing archive..."
( cd "$staging_dir" && sha256sum "$(basename "$RUN_TAR")" ) > "$RUN_TAR_SHA"

archive_bytes="$(du -b "$RUN_TAR" | awk '{print $1}')"
archive_sha256="$(awk '{print $1}' "$RUN_TAR_SHA")"

cat > "$RUN_META" <<EOF
output_dir: ${out_dir}
created_at: $(date -u +%Y-%m-%dT%H:%M:%SZ)
hostname: $(hostname)
user: $(id -un)
source_bytes: ${source_bytes}
archive_path: ${RUN_TAR}
archive_bytes: ${archive_bytes}
archive_sha256: ${archive_sha256}
manifest: $(basename "$RUN_SHA")
tar_version: $(tar --version 2>&1 | head -n 1)
zstd_version: $(zstd --version 2>&1 | head -n 1)
EOF

cat <<EOF

---------------------------------------------------------------------------
Staging complete.

  archive:   ${RUN_TAR}
  manifest:  ${RUN_SHA}
  metadata:  ${RUN_META}
  sha256:    ${archive_sha256}

Next steps (run from your LOCAL machine; edit <user>@<host>:<dest>):

  # 1. Pull the four sidecar files
  rsync -avP --partial \\
      <user>@<host>:${RUN_TAR} \\
      <user>@<host>:${RUN_TAR_SHA} \\
      <user>@<host>:${RUN_SHA} \\
      <user>@<host>:${RUN_META} \\
      <dest>/

  # 2. Verify the archive hash
  cd <dest> && sha256sum -c $(basename "$RUN_TAR_SHA")

  # 3. Extract (macOS: brew install zstd)
  tar -I zstd -xvf $(basename "$RUN_TAR")

  # 4. Verify every extracted file against the per-file manifest
  sha256sum -c $(basename "$RUN_SHA")
---------------------------------------------------------------------------
EOF
