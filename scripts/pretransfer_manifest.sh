#!/usr/bin/env bash
# pretransfer_manifest.sh — stage a NemaScan-sims-nf pipeline output directory
# as a reproducible tar.zst archive plus a per-file SHA-256 manifest and a
# human-readable run metadata file.
#
# Usage:
#   pretransfer_manifest.sh <output_dir> [staging_dir] [--force]
#
#   <output_dir>   Pipeline output directory (must contain a db/ subdir).
#                  Example: /vast/eande106/$USER/Sims_ce-cb-ct_nqtl1-50_h2grid_50reps
#   [staging_dir]  Where to write archive + sidecar files.
#                  Default: parent directory of <output_dir>.
#   --force        Re-create staging artifacts even if they are newer than the
#                  source. Without --force, the script is idempotent.
#
# Exit codes:
#   0  success
#   1  validation / usage error
#   2  staging directory write failure
#   3  tar / zstd failure
#   4  hash failure
#
# Primary target: Rockfish HPC (Linux, GNU coreutils + GNU tar >= 1.28 + zstd).
# Also runs on macOS for local smoke tests — BSD stat/du are detected and
# handled via small wrappers. macOS needs GNU tar available as `gtar`
# (`brew install gnu-tar`). The script hard-requires GNU tar because
# reproducible archive hashes depend on `--sort=name --mtime=...
# --owner=0 --group=0 --numeric-owner`, none of which BSD tar supports.

set -euo pipefail

err_trap() {
    local ec=$? line=${BASH_LINENO[0]:-?}
    printf 'pretransfer_manifest.sh: FAILED at line %s (exit %s)\n' "$line" "$ec" >&2
    exit "$ec"
}
trap err_trap ERR

usage() {
    sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

# ---- portable helpers ------------------------------------------------------
if stat --version >/dev/null 2>&1; then
    stat_size()  { stat -c%s "$1"; }
    stat_mtime() { stat -c%Y "$1"; }
else
    stat_size()  { stat -f%z "$1"; }
    stat_mtime() { stat -f%m "$1"; }
fi

dir_bytes() {
    # Portable: sum file sizes via find + stat_size.
    find "$1" -type f -print0 \
        | xargs -0 -n 50 bash -c 'for f; do stat_size "$f"; done' _ \
        | awk '{ s += $1 } END { print s+0 }'
}

newest_mtime_in() {
    find "$1" -type f -print0 \
        | xargs -0 -n 50 bash -c 'for f; do stat_mtime "$f"; done' _ \
        | sort -n | tail -1
}

# Require GNU tar (gtar on macOS-with-brew, tar on Linux). BSD tar lacks the
# flags we depend on for reproducibility, so fall back is NOT safe.
detect_gnu_tar() {
    if command -v gtar >/dev/null 2>&1 && gtar --version 2>&1 | head -1 | grep -q '^tar (GNU tar)'; then
        TAR=gtar
        return 0
    fi
    if command -v tar >/dev/null 2>&1 && tar --version 2>&1 | head -1 | grep -q '^tar (GNU tar)'; then
        TAR=tar
        return 0
    fi
    return 1
}
if ! detect_gnu_tar; then
    cat >&2 <<'NOTAR'
ERR: GNU tar (>= 1.28) is required for reproducible archives but was not found.
     On macOS:  brew install gnu-tar  (provides `gtar`)
     On Linux:  it should already be your default `tar` — check `tar --version`.
NOTAR
    exit 1
fi

export -f stat_size stat_mtime

# ---- argument parsing ------------------------------------------------------
force=0
positional=()
for arg in "$@"; do
    case "$arg" in
        --force) force=1 ;;
        -h|--help) usage ;;
        --) shift ;;
        *) positional+=("$arg") ;;
    esac
done

if [[ ${#positional[@]} -lt 1 || ${#positional[@]} -gt 2 ]]; then
    usage
fi

out_dir="${positional[0]%/}"
staging_dir="${positional[1]:-$(dirname "$out_dir")}"
staging_dir="${staging_dir%/}"

# ---- validation ------------------------------------------------------------
if [[ ! -d "$out_dir" ]]; then
    printf 'ERR: output dir %q does not exist or is not a directory\n' "$out_dir" >&2
    exit 1
fi
if [[ ! -d "$out_dir/db" ]]; then
    printf 'ERR: %q does not look like a pipeline output (missing db/ subdir)\n' "$out_dir" >&2
    exit 1
fi
if ! mkdir -p "$staging_dir"; then
    printf 'ERR: cannot create staging dir %q\n' "$staging_dir" >&2
    exit 2
fi
if [[ ! -w "$staging_dir" ]]; then
    printf 'ERR: staging dir %q is not writable\n' "$staging_dir" >&2
    exit 2
fi
for tool in "$TAR" zstd sha256sum find sort du; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        printf 'ERR: required tool %q not found in PATH\n' "$tool" >&2
        exit 1
    fi
done

base=$(basename "$out_dir")
parent=$(cd "$(dirname "$out_dir")" && pwd)
archive="$staging_dir/${base}.tar.zst"
archive_sha="$staging_dir/${base}.tar.zst.sha256"
manifest="$staging_dir/${base}.sha256"
meta="$staging_dir/${base}.manifest.txt"

# ---- idempotency guard -----------------------------------------------------
all_exist=1
for f in "$archive" "$archive_sha" "$manifest" "$meta"; do
    [[ -f "$f" ]] || { all_exist=0; break; }
done

if [[ $all_exist -eq 1 && $force -eq 0 ]]; then
    newest_source=$(newest_mtime_in "$out_dir")
    oldest_staging=$(
        for f in "$archive" "$archive_sha" "$manifest" "$meta"; do
            stat_mtime "$f"
        done | sort -n | head -1
    )
    if awk -v s="${newest_source:-0}" -v t="${oldest_staging:-0}" 'BEGIN { exit !(s <= t) }'; then
        printf 'pretransfer_manifest.sh: staging is current for %s\n' "$base"
        printf '  (use --force to regenerate)\n'
        printf '  archive: %s\n' "$archive"
        exit 0
    fi
fi

# ---- per-file SHA-256 manifest --------------------------------------------
printf 'Generating per-file manifest...\n'
(
    cd "$parent"
    find "$base" -type f -print0 \
        | LC_ALL=C sort -z \
        | xargs -0 sha256sum
) > "$manifest" || { printf 'ERR: manifest generation failed\n' >&2; exit 4; }

file_count=$(wc -l < "$manifest" | tr -d ' ')
source_bytes=$(dir_bytes "$out_dir")

# ---- reproducible tar.zst archive -----------------------------------------
printf 'Creating reproducible tar.zst archive (this may take several minutes)...\n'
(
    cd "$parent"
    "$TAR" --sort=name \
           --owner=0 --group=0 --numeric-owner \
           --mtime='UTC 2020-01-01' \
           -cf - "$base"
) | zstd -T0 -9 -q -o "$archive" \
    || { printf 'ERR: tar|zstd failed\n' >&2; exit 3; }

# ---- hash the archive ------------------------------------------------------
printf 'Hashing archive...\n'
(
    cd "$staging_dir"
    sha256sum "$(basename "$archive")"
) > "$archive_sha" || { printf 'ERR: archive hash failed\n' >&2; exit 4; }

archive_hash=$(awk '{print $1}' "$archive_sha")
archive_bytes=$(stat_size "$archive")

# ---- summary metadata ------------------------------------------------------
ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
script_dir=$(cd "$(dirname "$0")" && pwd)
pipeline_commit=$(git -C "$script_dir" rev-parse HEAD 2>/dev/null || echo "unknown")
tar_version=$("$TAR" --version 2>&1 | head -1)
zstd_version=$(zstd --version 2>&1 | head -1)
host=$(hostname)
whoami_val=$(id -un)

cat > "$meta" <<EOF
# pretransfer_manifest.sh — staging metadata
output_dir:       $out_dir
created_at:       $ts
hostname:         $host
user:             $whoami_val
pipeline_commit:  $pipeline_commit
file_count:       $file_count
source_bytes:     $source_bytes
archive_path:     $archive
archive_bytes:    $archive_bytes
archive_sha256:   $archive_hash
manifest:         $manifest
tar_version:      $tar_version
zstd_version:     $zstd_version
EOF

# ---- next-step instructions ------------------------------------------------
cat <<EOF

---------------------------------------------------------------------------
Staging complete.

  archive:   $archive  ($archive_bytes bytes)
  manifest:  $manifest  ($file_count files)
  metadata:  $meta
  sha256:    $archive_hash

Next steps (run from your LOCAL machine):

  # 1. Pull the four sidecar files
  rsync -avP --partial \\
    ${whoami_val}@${host}:${archive} \\
    ${whoami_val}@${host}:${archive_sha} \\
    ${whoami_val}@${host}:${manifest} \\
    ${whoami_val}@${host}:${meta} \\
    ./

  # 2. Verify the archive integrity (single hash, fast)
  sha256sum -c $(basename "$archive_sha")

  # 3. Extract (macOS: brew install zstd)
  tar -I zstd -xvf $(basename "$archive")

  # 4. Verify every extracted file against the per-file manifest
  sha256sum -c $(basename "$manifest")

Only after step 4 reports OK for every file is it safe to clean up scratch.
---------------------------------------------------------------------------
EOF
