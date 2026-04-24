#!/usr/bin/env bash
# scan_workdir.sh — Audit Nextflow work directory sizes on HPC
#
# Usage:
#   bash scan_workdir.sh <work_dir> <output_tsv> [trace_file]
#
# Arguments:
#   work_dir    Path to the Nextflow work directory
#               e.g. /scratch4/eande106/Ryan_tmp/nf-work-panelsize
#   output_tsv  Path for the output TSV file
#   trace_file  (optional) Path to a Nextflow execution trace file
#               e.g. Sims_panelsize_ce_.../pipeline_info/execution_trace_*.txt
#               If omitted, process names are inferred from .command.run headers.
#
# Output TSV columns:
#   work_hash           relative path to task dir (e.g. ab/cdef1234...)
#   process             Nextflow process name (e.g. PLINK_RECODE_MS_VCF)
#   tag                 task tag / meta.id (e.g. ce.n300.r29_0.05)
#   panel_size          integer panel size (100/200/300/400/500) or NA
#   total_bytes         total disk usage of the task dir in bytes
#   n_files             total file entries (real files + symlinks)
#   n_symlinks          number of symbolic links
#   n_real_files        number of regular (non-symlink) files
#   largest_file_bytes  size in bytes of the largest non-symlink file
#   largest_file_name   filename of the largest non-symlink file
#
# Performance:
#   Uses batched find+du and a single recursive find -printf inventory pass
#   to avoid 58K+ subprocess calls. Scales to 100K+ task directories.

set -euo pipefail

WORK_DIR="${1:?ERROR: work_dir argument required. Usage: $0 <work_dir> <output_tsv> [trace_file]}"
OUT_TSV="${2:?ERROR: output_tsv argument required. Usage: $0 <work_dir> <output_tsv> [trace_file]}"
TRACE_FILE="${3:-}"

WORK_DIR="${WORK_DIR%/}"   # strip trailing slash

echo "[scan] ──────────────────────────────────────────────────────"
echo "[scan] Nextflow work directory audit"
echo "[scan] Work dir:  $WORK_DIR"
echo "[scan] Output:    $OUT_TSV"
echo "[scan] Trace:     ${TRACE_FILE:-not provided}"
echo "[scan] ──────────────────────────────────────────────────────"

# ─── Temp files ──────────────────────────────────────────────────────────────
TMP_TRACE=$(mktemp /tmp/nf_scan_trace.XXXXXX)
TMP_DU=$(mktemp /tmp/nf_scan_du.XXXXXX)
TMP_INV=$(mktemp /tmp/nf_scan_inv.XXXXXX)
TMP_METRICS=$(mktemp /tmp/nf_scan_metrics.XXXXXX)

cleanup() { rm -f "$TMP_TRACE" "$TMP_DU" "$TMP_INV" "$TMP_METRICS"; }
trap cleanup EXIT

# ─── STEP 1: Parse trace file ────────────────────────────────────────────────
# Trace TSV columns (header row present):
#   task_id  hash  native_id  name  status  exit  submit  ...
# hash format: "ab/cdef12"  (2-char prefix + first 6 chars of full hash)
# name format: "[WORKFLOW:]PROCESS_NAME [(tag)]"

echo "[scan] [1/5] Parsing trace file..."
if [[ -n "$TRACE_FILE" && -f "$TRACE_FILE" ]]; then
    awk -F'\t' '
    NR == 1 { next }           # skip header
    NF < 4  { next }           # skip malformed rows
    {
        hash = $2
        name = $4

        # Strip workflow prefix: "WORKFLOW:PROCESS" -> "PROCESS"
        sub(/^.*:/, "", name)

        # Extract tag from parentheses, then strip from name
        tag = ""
        if (match(name, /\(([^)]+)\)/)) {
            tag = substr(name, RSTART+1, RLENGTH-2)
            sub(/ *\([^)]+\)$/, "", name)
        }

        # Trim whitespace
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", name)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", tag)

        if (name != "") print hash "\t" name "\t" tag
    }' "$TRACE_FILE" > "$TMP_TRACE"
    echo "[scan]   Loaded $(wc -l < "$TMP_TRACE") trace entries"
else
    touch "$TMP_TRACE"
    echo "[scan]   No trace file — will use .command.run fallback for all dirs"
fi

# ─── STEP 2: Batch du ─────────────────────────────────────────────────────────
# Single xargs-batched du call; outputs: rel_key TAB total_bytes
echo "[scan] [2/5] Batch du across all task dirs (this may take several minutes)..."
(
    cd "$WORK_DIR"
    find . -mindepth 2 -maxdepth 2 -type d -print0 \
        | xargs -0 du --block-size=1 --summarize 2>/dev/null \
        | awk '{ sub(/^\.\//, "", $2); print $2 "\t" $1 }'
) > "$TMP_DU"
NTASKS=$(wc -l < "$TMP_DU")
echo "[scan]   du complete: $NTASKS task dirs found"

# ─── STEP 3: Single find pass for file inventory ─────────────────────────────
# Emit: type(f/l) TAB size_bytes TAB path_relative_to_WORK_DIR
# -mindepth 3 skips the 2-char and hash dir entries themselves
echo "[scan] [3/5] Collecting file inventory (single find pass)..."
find "$WORK_DIR" -mindepth 3 \
    -not -path "${WORK_DIR}/.nextflow*" \
    \( -type f -o -type l \) \
    -printf '%y\t%s\t%P\n' \
> "$TMP_INV"
echo "[scan]   Inventory: $(wc -l < "$TMP_INV") files+symlinks"

# ─── STEP 4: Aggregate per-task-dir metrics ───────────────────────────────────
# Groups inventory by the first two path components (the task dir rel key).
# Accumulates: n_files, n_symlinks, n_real_files, largest_file_bytes, largest_file_name
echo "[scan] [4/5] Aggregating per-dir metrics..."
awk -F'\t' '
{
    path = $3
    n = split(path, parts, "/")
    if (n < 2) next
    key   = parts[1] "/" parts[2]
    type  = $1
    size  = $2 + 0
    fname = parts[n]

    n_files[key]++
    if (type == "l") {
        n_sym[key]++
    } else {
        n_real[key]++
        if (size > max_sz[key]) {
            max_sz[key] = size
            max_fn[key] = fname
        }
    }
}
END {
    for (k in n_files) {
        ns = (k in n_sym)  ? n_sym[k]  : 0
        nr = (k in n_real) ? n_real[k] : 0
        ms = (k in max_sz) ? max_sz[k] : 0
        mn = (k in max_fn) ? max_fn[k] : ""
        print k "\t" n_files[k] "\t" ns "\t" nr "\t" ms "\t" mn
    }
}' "$TMP_INV" > "$TMP_METRICS"
echo "[scan]   Metrics aggregated: $(wc -l < "$TMP_METRICS") dirs"

# ─── STEP 5: Join all sources and write output TSV ────────────────────────────
# Joins: trace (short_hash -> process, tag)
#        du    (rel_key    -> total_bytes)
#        metrics (rel_key  -> n_files, n_sym, n_real, max_bytes, max_name)
#
# Panel size extracted from tag via pattern .n{DIGITS}[._]
# (handles: ce.n300.r29_0.05, ce.n100.r01, etc.)
echo "[scan] [5/5] Joining data and writing output TSV..."

printf 'work_hash\tprocess\ttag\tpanel_size\ttotal_bytes\tn_files\tn_symlinks\tn_real_files\tlargest_file_bytes\tlargest_file_name\n' > "$OUT_TSV"

awk -F'\t' \
    -v trace_file="$TMP_TRACE" \
    -v du_file="$TMP_DU" \
'
# File 1: trace  (short_hash TAB process TAB tag)
FILENAME == trace_file {
    if (NF >= 3) {
        trace_proc[$1] = $2
        trace_tag[$1]  = $3
    }
    next
}

# File 2: du map  (rel_key TAB total_bytes)
FILENAME == du_file {
    du[$1] = $2
    next
}

# File 3: metrics  (rel_key TAB n_files TAB n_sym TAB n_real TAB max_bytes TAB max_name)
{
    rel_key   = $1
    n_files   = $2
    n_sym     = $3
    n_real    = $4
    max_bytes = $5
    max_fn    = $6

    # Build short hash from rel_key: "ab" + "/" + first 6 chars of second component
    n = split(rel_key, parts, "/")
    short_hash = (n >= 2) ? (parts[1] "/" substr(parts[2], 1, 6)) : ""

    process     = (short_hash in trace_proc) ? trace_proc[short_hash] : ""
    tag         = (short_hash in trace_tag)  ? trace_tag[short_hash]  : ""
    total_bytes = (rel_key in du) ? du[rel_key] : 0

    # Extract panel size from tag: match .n{digits}[._\0]
    panel_size = "NA"
    if (tag != "") {
        tmp = tag
        if (sub(/.*\.n/, "", tmp) && match(tmp, /^[0-9]+/)) {
            panel_size = substr(tmp, 1, RLENGTH)
        }
    }
    # Fallback: try from largest filename
    if (panel_size == "NA" && max_fn != "") {
        tmp = max_fn
        if (sub(/.*\.n/, "", tmp) && match(tmp, /^[0-9]+/)) {
            panel_size = substr(tmp, 1, RLENGTH)
        }
    }

    print rel_key "\t" process "\t" tag "\t" panel_size "\t" \
          total_bytes "\t" n_files "\t" n_sym "\t" n_real "\t" \
          max_bytes "\t" max_fn
}
' "$TMP_TRACE" "$TMP_DU" "$TMP_METRICS" >> "$OUT_TSV"

# ─── Fallback: resolve empty process names via .command.run ──────────────────
# Only fires for LOCAL tasks or dirs missing from the trace (typically ≤300 rows).
MISSING=$(awk -F'\t' 'NR > 1 && $2 == ""' "$OUT_TSV" | wc -l)
if [[ "$MISSING" -gt 0 ]]; then
    echo "[scan] Resolving $MISSING dirs via .command.run fallback..."
    TMP_RESOLVED=$(mktemp /tmp/nf_scan_resolved.XXXXXX)
    head -1 "$OUT_TSV" > "$TMP_RESOLVED"

    while IFS=$'\t' read -r work_hash process tag panel_size total_bytes \
                               n_files n_sym n_real max_bytes max_fn; do
        if [[ -z "$process" ]]; then
            cmd_run="$WORK_DIR/$work_hash/.command.run"
            if [[ -f "$cmd_run" ]]; then
                process=$(grep -m1 '#SBATCH --job-name' "$cmd_run" 2>/dev/null \
                          | grep -oP 'nf-\K[A-Z_0-9]+' || true)
            fi
            [[ -z "$process" ]] && process="LOCAL"

            # Try to find tag from non-hidden filenames if still empty
            if [[ -z "$tag" ]]; then
                tag=$(ls "$WORK_DIR/$work_hash" 2>/dev/null \
                      | grep -v '^\.' \
                      | grep -oP '(?:ce|cb|ct)\.n[0-9]+\.[^[:space:]]+' \
                      | head -1 || true)
                if [[ -n "$tag" ]]; then
                    panel_size=$(printf '%s' "$tag" | grep -oP '(?<=\.n)\d+' | head -1 || echo "NA")
                fi
            fi
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$work_hash" "$process" "$tag" "$panel_size" \
            "$total_bytes" "$n_files" "$n_sym" "$n_real" \
            "$max_bytes" "$max_fn"
    done < <(tail -n +2 "$OUT_TSV") >> "$TMP_RESOLVED"

    mv "$TMP_RESOLVED" "$OUT_TSV"
    echo "[scan]   Fallback resolution complete"
fi

TOTAL=$(( $(wc -l < "$OUT_TSV") - 1 ))
RESOLVED=$(awk -F'\t' 'NR > 1 && $2 != "" && $2 != "LOCAL"' "$OUT_TSV" | wc -l)
UNRESOLVED=$(awk -F'\t' 'NR > 1 && ($2 == "" || $2 == "LOCAL")' "$OUT_TSV" | wc -l)

echo "[scan] ──────────────────────────────────────────────────────"
echo "[scan] Complete."
echo "[scan]   Total task dirs:  $TOTAL"
echo "[scan]   Named processes:  $RESOLVED"
echo "[scan]   LOCAL/unresolved: $UNRESOLVED"
echo "[scan]   Output:           $OUT_TSV"
echo "[scan] ──────────────────────────────────────────────────────"
