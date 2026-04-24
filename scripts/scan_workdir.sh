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
#               If omitted, process names are extracted from .command.run files
#               via a single batched grep pass (NOT a per-dir shell loop).
#
# Output TSV columns:
#   work_hash           relative path to task dir (e.g. ab/cdef1234...)
#   process             Nextflow process name (e.g. PLINK_RECODE_MS_VCF)
#   tag                 task tag / meta.id (e.g. ce.n300.r29_0.05)
#   panel_size          integer panel size extracted from tag, or NA
#   total_bytes         total disk usage of the task dir in bytes
#   n_files             total file entries (real files + symlinks)
#   n_symlinks          number of symbolic links
#   n_real_files        number of regular (non-symlink) files
#   largest_file_bytes  size in bytes of the largest non-symlink file
#   largest_file_name   filename of the largest non-symlink file

set -euo pipefail

WORK_DIR="${1:?ERROR: work_dir argument required. Usage: $0 <work_dir> <output_tsv> [trace_file]}"
OUT_TSV="${2:?ERROR: output_tsv argument required. Usage: $0 <work_dir> <output_tsv> [trace_file]}"
TRACE_FILE="${3:-}"

WORK_DIR="${WORK_DIR%/}"

echo "[scan] ──────────────────────────────────────────────────────"
echo "[scan] Nextflow work directory audit"
echo "[scan] Work dir:  $WORK_DIR"
echo "[scan] Output:    $OUT_TSV"
echo "[scan] Trace:     ${TRACE_FILE:-not provided}"
echo "[scan] Started:   $(date)"
echo "[scan] ──────────────────────────────────────────────────────"

# ─── Temp files ──────────────────────────────────────────────────────────────
TMP_TRACE=$(mktemp /tmp/nf_scan_trace.XXXXXX)
TMP_DU=$(mktemp /tmp/nf_scan_du.XXXXXX)
TMP_INV=$(mktemp /tmp/nf_scan_inv.XXXXXX)
TMP_METRICS=$(mktemp /tmp/nf_scan_metrics.XXXXXX)
TMP_CMD_MAP=$(mktemp /tmp/nf_scan_cmdmap.XXXXXX)

cleanup() { rm -f "$TMP_TRACE" "$TMP_DU" "$TMP_INV" "$TMP_METRICS" "$TMP_CMD_MAP"; }
trap cleanup EXIT

# ─── Helper: run a command in background and print progress every 30s ─────────
run_with_progress() {
    local label="$1"; shift
    local outfile="$1"; shift
    # Remaining args are the command to run
    "$@" > "$outfile" &
    local pid=$!
    local elapsed=0
    while kill -0 "$pid" 2>/dev/null; do
        sleep 30
        elapsed=$(( elapsed + 30 ))
        local lines
        lines=$(wc -l < "$outfile" 2>/dev/null || echo "?")
        echo "[scan]   ... $label: ${elapsed}s elapsed, ~${lines} lines written"
    done
    wait "$pid"
}

# ─── STEP 1: Parse trace file ────────────────────────────────────────────────
echo "[scan] [1/5] Parsing trace file..."
if [[ -n "$TRACE_FILE" && -f "$TRACE_FILE" ]]; then
    awk -F'\t' '
    NR == 1 { next }
    NF < 4  { next }
    {
        hash = $2
        name = $4
        sub(/^.*:/, "", name)          # strip workflow prefix
        tag = ""
        if (match(name, /\(([^)]+)\)/)) {
            tag = substr(name, RSTART+1, RLENGTH-2)
            sub(/ *\([^)]+\)$/, "", name)
        }
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", name)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", tag)
        if (name != "") print hash "\t" name "\t" tag
    }' "$TRACE_FILE" > "$TMP_TRACE"
    echo "[scan]   Loaded $(wc -l < "$TMP_TRACE") trace entries"
else
    touch "$TMP_TRACE"
    echo "[scan]   No trace file — will extract process names from .command.run (step 5)"
fi

# ─── STEP 2: Batch du ─────────────────────────────────────────────────────────
# Outputs: rel_key TAB total_bytes
# Runs in background with 30s progress ticks.
echo "[scan] [2/5] Batch du across all task dirs..."
echo "[scan]   (GPFS/Lustre can take 10–40 min for 58K dirs — progress every 30s)"

run_with_progress "dirs in du" "$TMP_DU" \
    bash -c "cd '$WORK_DIR' && \
             find . -mindepth 2 -maxdepth 2 -type d -print0 \
             | xargs -0 du --block-size=1 --summarize 2>/dev/null \
             | awk '{sub(/^\\.\\//,\"\",\$2); print \$2\"\\t\"\$1}'"

NTASKS=$(wc -l < "$TMP_DU")
echo "[scan]   du complete: $NTASKS task dirs found"

# ─── STEP 3: Single find pass for file inventory ─────────────────────────────
# Outputs: type(f/l) TAB size_bytes TAB path_relative_to_WORK_DIR
echo "[scan] [3/5] Collecting file inventory (single find pass)..."
echo "[scan]   (stat-ing all files in $NTASKS dirs — progress every 30s)"

run_with_progress "files inventoried" "$TMP_INV" \
    find "$WORK_DIR" -mindepth 3 \
        -not -path "${WORK_DIR}/.nextflow*" \
        \( -type f -o -type l \) \
        -printf '%y\t%s\t%P\n'

echo "[scan]   Inventory: $(wc -l < "$TMP_INV") files+symlinks"

# ─── STEP 4: Aggregate per-task-dir metrics ───────────────────────────────────
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

# ─── STEP 5: Build process map ────────────────────────────────────────────────
# Primary source: trace file (already loaded into TMP_TRACE).
# Fallback: single batched grep over all .command.run files — NOT a shell loop.
# Produces TMP_CMD_MAP: rel_key TAB process  (for SLURM tasks, via SBATCH job-name)
echo "[scan] [5/5] Building process map and writing output TSV..."

WDIR_LEN=${#WORK_DIR}

if [[ -s "$TMP_TRACE" ]]; then
    echo "[scan]   Process names from trace file (no .command.run scan needed)"
    touch "$TMP_CMD_MAP"
else
    echo "[scan]   No trace — scanning all .command.run files (batched grep, one pass)..."
    find "$WORK_DIR" -mindepth 3 -maxdepth 3 -name '.command.run' -print0 \
        | xargs -0 grep -H '#SBATCH --job-name' 2>/dev/null \
        | awk -v wdir_len="$WDIR_LEN" '{
            # Input: /full/path/.command.run:#SBATCH --job-name nf-PROCESSNAME
            colon = index($0, ":")
            filepath = substr($0, 1, colon - 1)
            content  = substr($0, colon + 1)
            # Strip work dir prefix (+2 to skip trailing slash) and /.command.run suffix
            rel_path = substr(filepath, wdir_len + 2)
            sub("/\\.command\\.run$", "", rel_path)
            # Extract process name (text after "nf-")
            if (match(content, /nf-[A-Z_0-9]+/)) {
                process = substr(content, RSTART + 3, RLENGTH - 3)
            } else {
                process = ""
            }
            if (rel_path != "" && process != "") print rel_path "\t" process
        }' > "$TMP_CMD_MAP"
    echo "[scan]   .command.run map: $(wc -l < "$TMP_CMD_MAP") SLURM entries found"
fi

# ─── Join all sources → output TSV ───────────────────────────────────────────
printf 'work_hash\tprocess\ttag\tpanel_size\ttotal_bytes\tn_files\tn_symlinks\tn_real_files\tlargest_file_bytes\tlargest_file_name\n' > "$OUT_TSV"

awk -F'\t' \
    -v trace_file="$TMP_TRACE" \
    -v du_file="$TMP_DU" \
    -v cmd_file="$TMP_CMD_MAP" \
'
# File 1: trace  (short_hash TAB process TAB tag)
FILENAME == trace_file {
    if (NF >= 3) { trace_proc[$1] = $2; trace_tag[$1] = $3 }
    next
}
# File 2: du map  (rel_key TAB total_bytes)
FILENAME == du_file {
    du[$1] = $2
    next
}
# File 3: .command.run map  (rel_key TAB process)
FILENAME == cmd_file {
    if (NF >= 2) cmd_proc[$1] = $2
    next
}
# File 4: metrics  (rel_key TAB n_files TAB n_sym TAB n_real TAB max_bytes TAB max_name)
{
    rel_key   = $1
    n_files   = $2
    n_sym     = $3
    n_real    = $4
    max_bytes = $5
    max_fn    = $6

    # Short hash for trace lookup: "ab" + "/" + first 6 chars of second component
    n = split(rel_key, parts, "/")
    short_hash = (n >= 2) ? (parts[1] "/" substr(parts[2], 1, 6)) : ""

    # Process: trace > .command.run map > LOCAL
    if (short_hash in trace_proc) {
        process = trace_proc[short_hash]
        tag     = trace_tag[short_hash]
    } else if (rel_key in cmd_proc) {
        process = cmd_proc[rel_key]
        tag     = ""
    } else {
        process = "LOCAL"
        tag     = ""
    }

    total_bytes = (rel_key in du) ? du[rel_key] : 0

    # Extract panel size from tag: .n{digits}[._] pattern
    panel_size = "NA"
    if (tag != "") {
        tmp = tag
        if (sub(/.*\.n/, "", tmp) && match(tmp, /^[0-9]+/)) {
            panel_size = substr(tmp, 1, RLENGTH)
        }
    }
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
' "$TMP_TRACE" "$TMP_DU" "$TMP_CMD_MAP" "$TMP_METRICS" >> "$OUT_TSV"

TOTAL=$(( $(wc -l < "$OUT_TSV") - 1 ))
NAMED=$(awk -F'\t' 'NR>1 && $2!="" && $2!="LOCAL"' "$OUT_TSV" | wc -l)
LOCAL=$(awk -F'\t' 'NR>1 && $2=="LOCAL"' "$OUT_TSV" | wc -l)
EMPTY=$(awk -F'\t' 'NR>1 && $2==""' "$OUT_TSV" | wc -l)

echo "[scan] ──────────────────────────────────────────────────────"
echo "[scan] Complete.  $(date)"
echo "[scan]   Total task dirs:  $TOTAL"
echo "[scan]   Named processes:  $NAMED"
echo "[scan]   LOCAL tasks:      $LOCAL"
echo "[scan]   Unresolved:       $EMPTY"
echo "[scan]   Output:           $OUT_TSV"
echo "[scan] ──────────────────────────────────────────────────────"
