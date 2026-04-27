#!/usr/bin/env bash
# scan_workdir.sh — Audit Nextflow work directory sizes on HPC
#
# Usage:
#   bash scan_workdir.sh <work_dir> <output_tsv> <run_name> [OPTIONS]
#
# Arguments:
#   work_dir    Path to the Nextflow work directory
#               e.g. /scratch4/eande106/Ryan_tmp/nf-work-panelsize
#   output_tsv  Path for the output TSV file
#   run_name    Nextflow run name (e.g. spontaneous_swirles).
#               Run 'nextflow log' from your pipeline directory to list names.
#
# Options:
#   --processes PROC1,PROC2,...   Comma-separated process names to scan.
#                                 If omitted, all processes are included.
#   --nf_dir PATH                 Directory containing .nextflow/ history
#                                 (default: current working directory).
#   --parallel N                  Parallel du/find workers (default: 4).
#
# Output TSV columns:
#   work_hash           relative path to task dir (e.g. ab/cdef1234...)
#   process             Nextflow process name (e.g. PLINK_RECODE_MS_VCF)
#   tag                 task tag / meta.id (e.g. ce.n300.r29_0.05)
#   panel_size          integer panel size extracted from tag, or NA
#   total_bytes         total disk usage of the task dir in bytes (du --block-size=1)
#   n_files             total file entries (real files + symlinks)
#   n_symlinks          number of symbolic links
#   n_real_files        number of regular (non-symlink) files
#   largest_file_bytes  size in bytes of the largest non-symlink file
#   largest_file_name   filename of the largest non-symlink file

set -euo pipefail

WORK_DIR="${1:?ERROR: work_dir argument required. Usage: $0 <work_dir> <output_tsv> <run_name> [OPTIONS]}"
OUT_TSV="${2:?ERROR: output_tsv argument required. Usage: $0 <work_dir> <output_tsv> <run_name> [OPTIONS]}"
RUN_NAME="${3:?ERROR: run_name argument required (e.g. spontaneous_swirles). Run 'nextflow log' to list names.}"
shift 3

PROC_FILTER=""
NF_DIR="."
PARALLEL=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --processes) PROC_FILTER="$2"; shift 2 ;;
        --nf_dir)    NF_DIR="$2";      shift 2 ;;
        --parallel)  PARALLEL="$2";    shift 2 ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
    esac
done

WORK_DIR="${WORK_DIR%/}"
NF_DIR="${NF_DIR%/}"

echo "[scan] ──────────────────────────────────────────────────────"
echo "[scan] Nextflow work directory audit"
echo "[scan] Work dir:  $WORK_DIR"
echo "[scan] Run name:  $RUN_NAME"
echo "[scan] Output:    $OUT_TSV"
echo "[scan] Processes: ${PROC_FILTER:-all}"
echo "[scan] NF dir:    $NF_DIR"
echo "[scan] Parallel:  $PARALLEL workers"
echo "[scan] Started:   $(date)"
echo "[scan] ──────────────────────────────────────────────────────"

# ─── Temp files ──────────────────────────────────────────────────────────────
TMP_TASKS=$(mktemp /tmp/nf_scan_tasks.XXXXXX)
TMP_DU=$(mktemp /tmp/nf_scan_du.XXXXXX)
TMP_INV=$(mktemp /tmp/nf_scan_inv.XXXXXX)

cleanup() { rm -f "$TMP_TASKS" "$TMP_DU" "$TMP_INV"; }
trap cleanup EXIT

# ─── Helper: run a command in background and print progress every 30s ─────────
run_with_progress() {
    local label="$1"; shift
    local outfile="$1"; shift
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

# ─── STEP 1: nextflow log → task list ────────────────────────────────────────
# Outputs: process TAB tag TAB workdir  (no header, one row per task)
echo "[scan] [1/4] Fetching task list from 'nextflow log $RUN_NAME'..."
(cd "$NF_DIR" && nextflow log "$RUN_NAME" -f process,tag,workdir) > "$TMP_TASKS"
NTOTAL=$(wc -l < "$TMP_TASKS")
echo "[scan]   $NTOTAL total tasks found"

if [[ -n "$PROC_FILTER" ]]; then
    awk -F'\t' -v procs="$PROC_FILTER" '
        BEGIN { n = split(procs, a, ","); for (i = 1; i <= n; i++) f[a[i]] = 1 }
        $1 in f
    ' "$TMP_TASKS" > "${TMP_TASKS}.filt" && mv "${TMP_TASKS}.filt" "$TMP_TASKS"
    NTASKS=$(wc -l < "$TMP_TASKS")
    echo "[scan]   Filtered to '$PROC_FILTER': $NTASKS tasks ($(( NTOTAL - NTASKS )) excluded)"
else
    NTASKS=$NTOTAL
fi

if [[ "$NTASKS" -eq 0 ]]; then
    echo "[scan] ERROR: No matching tasks. Check run name and --processes filter." >&2
    exit 1
fi

# ─── STEP 2: du on task work dirs ─────────────────────────────────────────────
# Outputs: total_bytes TAB workdir  (du --block-size=1 format)
echo "[scan] [2/4] Measuring disk usage for $NTASKS dirs..."
echo "[scan]   (GPFS/Lustre can be slow — progress every 30s; $PARALLEL parallel workers)"

run_with_progress "dirs measured" "$TMP_DU" \
    bash -c "awk -F'\t' '{print \$3}' '$TMP_TASKS' \
        | xargs -P$PARALLEL du --block-size=1 --summarize 2>/dev/null"

echo "[scan]   du complete: $(wc -l < "$TMP_DU") dirs measured"

# ─── STEP 3: find file inventory ──────────────────────────────────────────────
# Outputs: type(f/l) TAB size_bytes TAB absolute_path
# Passes dirs in batches of 50 to a single find per invocation (-n50),
# running N batches in parallel. 'sh -c "find "$@"..." sh' puts dirs
# before find options so xargs appends them correctly.
echo "[scan] [3/4] Collecting file inventory for $NTASKS dirs..."
echo "[scan]   (progress every 30s)"

run_with_progress "files inventoried" "$TMP_INV" \
    bash -c "awk -F'\t' '{print \$3}' '$TMP_TASKS' \
        | xargs -P$PARALLEL -n50 sh -c \
            'find \"\$@\" -mindepth 1 \( -type f -o -type l \) -printf \"%y\t%s\t%p\n\"' sh"

echo "[scan]   Inventory: $(wc -l < "$TMP_INV") files+symlinks"

# ─── STEP 4: aggregate metrics + join all sources + write TSV ─────────────────
echo "[scan] [4/4] Aggregating metrics and writing output TSV..."

WDIR_LEN=${#WORK_DIR}

awk -F'\t' \
    -v tasks_file="$TMP_TASKS" \
    -v du_file="$TMP_DU" \
    -v wdir="$WORK_DIR" \
    -v wdir_len="$WDIR_LEN" \
'
# File 1: task list  (process TAB tag TAB workdir)
FILENAME == tasks_file {
    task_proc[$3] = $1
    task_tag[$3]  = $2
    next
}
# File 2: du output  (total_bytes TAB workdir)
FILENAME == du_file {
    du_bytes[$2] = $1 + 0
    next
}
# File 3: file inventory  (type TAB size TAB absolute_path)
{
    ftype = $1
    fsize = $2 + 0
    fpath = $3
    fname = fpath; sub(/.*\//, "", fname)

    # Derive absolute workdir: strip WORK_DIR prefix and keep first 2 path components
    rel = substr(fpath, wdir_len + 2)
    n = split(rel, parts, "/")
    if (n < 2) next
    wabs = wdir "/" parts[1] "/" parts[2]

    n_files[wabs]++
    if (ftype == "l") {
        n_sym[wabs]++
    } else {
        n_real[wabs]++
        if (fsize > max_sz[wabs]) {
            max_sz[wabs] = fsize
            max_fn[wabs] = fname
        }
    }
}
END {
    print "work_hash\tprocess\ttag\tpanel_size\ttotal_bytes\tn_files\tn_symlinks\tn_real_files\tlargest_file_bytes\tlargest_file_name"
    for (wabs in n_files) {
        process = (wabs in task_proc) ? task_proc[wabs] : "UNKNOWN"
        tag     = (wabs in task_tag)  ? task_tag[wabs]  : ""
        tb      = (wabs in du_bytes)  ? du_bytes[wabs]  : 0
        ns      = (wabs in n_sym)     ? n_sym[wabs]     : 0
        nr      = (wabs in n_real)    ? n_real[wabs]    : 0
        ms      = (wabs in max_sz)    ? max_sz[wabs]    : 0
        mn      = (wabs in max_fn)    ? max_fn[wabs]    : ""

        rel_key = substr(wabs, wdir_len + 2)

        # Extract panel size from tag: .n{digits}[._] pattern
        panel_size = "NA"
        tmp = tag
        if (sub(/.*\.n/, "", tmp) && match(tmp, /^[0-9]+/)) {
            panel_size = substr(tmp, 1, RLENGTH)
        }
        if (panel_size == "NA" && mn != "") {
            tmp = mn
            if (sub(/.*\.n/, "", tmp) && match(tmp, /^[0-9]+/)) {
                panel_size = substr(tmp, 1, RLENGTH)
            }
        }

        print rel_key "\t" process "\t" tag "\t" panel_size "\t" \
              tb "\t" n_files[wabs] "\t" ns "\t" nr "\t" ms "\t" mn
    }
}
' "$TMP_TASKS" "$TMP_DU" "$TMP_INV" > "$OUT_TSV"

TOTAL=$(( $(wc -l < "$OUT_TSV") - 1 ))
NAMED=$(awk -F'\t' 'NR>1 && $2 != "UNKNOWN" && $2 != ""' "$OUT_TSV" | wc -l)
UNKNOWN=$(awk -F'\t' 'NR>1 && ($2 == "UNKNOWN" || $2 == "")' "$OUT_TSV" | wc -l)

echo "[scan] ──────────────────────────────────────────────────────"
echo "[scan] Complete.  $(date)"
echo "[scan]   Total task dirs:  $TOTAL"
echo "[scan]   Named processes:  $NAMED"
echo "[scan]   Unresolved:       $UNKNOWN"
echo "[scan]   Output:           $OUT_TSV"
echo "[scan] ──────────────────────────────────────────────────────"
