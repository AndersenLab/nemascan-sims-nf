#!/usr/bin/env bash
#
# nemascan_run_with_cme_retry.sh
#
# Wrapper around `nextflow run` that auto-resumes the pipeline when it dies
# from the known ConcurrentModificationException at
# `BashWrapperBuilder.createContainerBuilder(BashWrapperBuilder.groovy:706)`
# (issue #175 — see
#  issues/175-gcta-concurrent-mod-exception/2026-05-19-rockfish-cme-post-cowal-rca.qmd
#  for the mechanism). The race is structural in Nextflow 24.10.1. The earlier
#  `CopyOnWriteArrayList` envWhitelist mitigation in `conf/rockfish.config` was
#  removed because it grew unbounded and OOM'd the driver at scale
#  (issues/new-unexpected-error/rca.md); this wrapper is now the primary defense,
#  detecting the abort and re-launching with `-resume`.
#
# All other failure modes (real bugs, SLURM time-limit kills, container
# pull failures, etc.) propagate to the caller WITHOUT retry. The detection
# is narrow on purpose so an unrelated bug cannot drive an infinite loop.
#
# Usage:
#   scripts/nemascan_run_with_cme_retry.sh [OPTIONS] -- <args passed to `nextflow run`>
#
# Options:
#   --max-retries N    Auto-resume attempts after the initial run. Default: 3.
#                      0 disables retries (the wrapper still runs once).
#   --nextflow CMD     Override the `nextflow` executable (default: `nextflow`).
#   --log-file PATH    Path to the .nextflow.log to scan after each attempt.
#                      Default: ./.nextflow.log (Nextflow's default location).
#   --no-resume-first  Do NOT pass `-resume` on the first attempt. By default
#                      `-resume` is added every attempt — harmless on a fresh
#                      run and required for continuing a prior aborted run.
#   -h, --help         Show this help and exit.
#
# Example (inside sbatch on Rockfish — see scripts/run_with_cme_retry.sbatch):
#   scripts/nemascan_run_with_cme_retry.sh --max-retries 3 -- \
#       main.nf -profile rockfish \
#       --strainfile strains.tsv --vcf 20231213 \
#       --output_dir /vast/eande106/${USER}/sim-$(date +%F)
#
# Where to run on Rockfish:
#   - Submit via `sbatch scripts/run_with_cme_retry.sbatch …` (a batch
#     allocation). Nextflow's slurm executor then submits child tasks
#     normally from inside that allocation.
#   - Do NOT run this wrapper from an interactive `srun --pty` shell:
#     Rockfish does not allow interactive srun allocations to submit
#     sbatch jobs, so Nextflow's slurm executor would fail to dispatch
#     child tasks. (Login node + sbatch submission is the supported
#     path.)
#
# Exit codes:
#   0        Pipeline exited 0 on some attempt.
#   non-zero Last Nextflow exit code, whether the abort was a detected CME
#            with retries exhausted or any other failure (no retry).

set -uo pipefail

# --------------------------------------------------------------------------
# Defaults & arg parsing
# --------------------------------------------------------------------------

MAX_RETRIES=3
NEXTFLOW_CMD="nextflow"
LOG_FILE=".nextflow.log"
RESUME_ON_FIRST=1

print_help() {
    sed -n '2,40p' "$0" | sed -e 's/^# \{0,1\}//'
}

# Pull wrapper options off the front; everything after `--` is for Nextflow.
nf_args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-retries)
            MAX_RETRIES="$2"; shift 2 ;;
        --max-retries=*)
            MAX_RETRIES="${1#*=}"; shift ;;
        --nextflow)
            NEXTFLOW_CMD="$2"; shift 2 ;;
        --nextflow=*)
            NEXTFLOW_CMD="${1#*=}"; shift ;;
        --log-file)
            LOG_FILE="$2"; shift 2 ;;
        --log-file=*)
            LOG_FILE="${1#*=}"; shift ;;
        --no-resume-first)
            RESUME_ON_FIRST=0; shift ;;
        -h|--help)
            print_help; exit 0 ;;
        --)
            shift; nf_args=("$@"); break ;;
        *)
            echo "ERR: unknown wrapper option '$1' (did you forget '--' before the nextflow args?)" >&2
            exit 2 ;;
    esac
done

if [[ ${#nf_args[@]} -eq 0 ]]; then
    echo "ERR: no nextflow arguments supplied. Use '-- main.nf -profile rockfish ...'" >&2
    print_help >&2
    exit 2
fi

if ! [[ "$MAX_RETRIES" =~ ^[0-9]+$ ]]; then
    echo "ERR: --max-retries must be a non-negative integer (got '$MAX_RETRIES')" >&2
    exit 2
fi

if ! command -v "$NEXTFLOW_CMD" >/dev/null 2>&1; then
    echo "ERR: '$NEXTFLOW_CMD' not on PATH. Activate the nextflow conda env first (e.g. 'conda activate /data/eande106/software/conda_envs/nf24_env')." >&2
    exit 2
fi

# --------------------------------------------------------------------------
# CME-at-706 detection
# --------------------------------------------------------------------------
#
# We match the deepest, most specific frame from the issue #175 stack:
#
#   java.util.ConcurrentModificationException: null
#       at java.base/java.util.ArrayList$Itr.checkForComodification(ArrayList.java:...)
#       at java.base/java.util.ArrayList$Itr.next(ArrayList.java:...)
#       at nextflow.executor.BashWrapperBuilder.createContainerBuilder(BashWrapperBuilder.groovy:706)
#
# Matching on the BashWrapperBuilder.groovy:706 substring alone is enough:
# the line number + class + method triple is unique across the NF 24.10.1
# codebase, and the log frame format is stable. We only look at the last
# 500 lines so that historical CMEs from a prior attempt (e.g. when the log
# wasn't rotated between runs) do not produce a false positive — the
# terminating frame is always near the end of the log.

CME_FRAME_NEEDLE='BashWrapperBuilder.createContainerBuilder(BashWrapperBuilder.groovy:706)'

cme_at_706_in_log() {
    local logfile="$1"
    if [[ ! -f "$logfile" ]]; then
        return 1
    fi
    tail -n 500 "$logfile" | grep -qF "$CME_FRAME_NEEDLE"
}

# --------------------------------------------------------------------------
# Signal handling: forward SIGTERM/SIGINT to the running nextflow child so
# SLURM scancel and Ctrl-C exit cleanly instead of leaving an orphan.
# --------------------------------------------------------------------------

current_nf_pid=""
forward_signal() {
    local sig="$1"
    if [[ -n "$current_nf_pid" ]] && kill -0 "$current_nf_pid" 2>/dev/null; then
        echo
        echo "[wrapper] forwarding $sig to nextflow pid $current_nf_pid"
        kill "-$sig" "$current_nf_pid" 2>/dev/null || true
    fi
}
trap 'forward_signal TERM' TERM
trap 'forward_signal INT'  INT

# --------------------------------------------------------------------------
# Banner
# --------------------------------------------------------------------------

print_banner() {
    cat <<EOF
================================================================================
  nemascan_run_with_cme_retry.sh — issue #175 auto-resume wrapper
--------------------------------------------------------------------------------
  nextflow:        $NEXTFLOW_CMD ($("$NEXTFLOW_CMD" -version 2>/dev/null | sed -n 's/.*version \([^ ]*\).*/\1/p' | head -1 || echo "version: unknown"))
  log file:        $LOG_FILE
  max retries:     $MAX_RETRIES
  cwd:             $PWD
  hostname:        $(hostname)
  SLURM job id:    ${SLURM_JOB_ID:-not in SLURM}
  trigger:         $CME_FRAME_NEEDLE
  passed to NF:    ${nf_args[*]}
================================================================================
EOF
}

print_attempt_header() {
    local attempt="$1" total_budget="$2" use_resume="$3"
    cat <<EOF

--------------------------------------------------------------------------------
  Attempt $attempt of $total_budget$( [[ "$use_resume" == "1" ]] && echo "  (with -resume)" )
  Started: $(date -Iseconds)
--------------------------------------------------------------------------------
EOF
}

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

start_epoch=$(date +%s)
total_budget=$((1 + MAX_RETRIES))

# Per-attempt records, parallel arrays
attempt_codes=()
attempt_reasons=()

print_banner

final_exit=0
for ((attempt=1; attempt <= total_budget; attempt++)); do
    # First attempt obeys --no-resume-first; subsequent attempts always resume.
    if [[ $attempt -eq 1 ]]; then
        use_resume=$RESUME_ON_FIRST
    else
        use_resume=1
    fi

    run_args=("run")
    if [[ "$use_resume" == "1" ]]; then
        run_args+=("-resume")
    fi
    run_args+=("${nf_args[@]}")

    print_attempt_header "$attempt" "$total_budget" "$use_resume"
    echo "[wrapper] command: $NEXTFLOW_CMD ${run_args[*]}"
    echo

    # Run nextflow in background so we can capture its PID for signal forwarding.
    "$NEXTFLOW_CMD" "${run_args[@]}" &
    current_nf_pid=$!
    wait "$current_nf_pid"
    rc=$?
    current_nf_pid=""

    attempt_codes+=("$rc")
    echo
    echo "[wrapper] attempt $attempt exited with code $rc at $(date -Iseconds)"

    if [[ $rc -eq 0 ]]; then
        attempt_reasons+=("success")
        final_exit=0
        break
    fi

    if cme_at_706_in_log "$LOG_FILE"; then
        attempt_reasons+=("CME at BashWrapperBuilder:706 — auto-resume eligible")
        echo "[wrapper] detected CME at line 706 in $LOG_FILE"
        if [[ $attempt -lt $total_budget ]]; then
            echo "[wrapper] will retry with -resume ($((total_budget - attempt)) attempts left)"
            continue
        else
            echo "[wrapper] retry budget exhausted ($MAX_RETRIES retries used)"
            final_exit=$rc
            break
        fi
    else
        attempt_reasons+=("non-CME failure — not retrying")
        echo "[wrapper] failure is NOT the CME-at-706 pattern; propagating exit code"
        final_exit=$rc
        break
    fi
done

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

end_epoch=$(date +%s)
elapsed=$((end_epoch - start_epoch))

cat <<EOF

================================================================================
  Summary
--------------------------------------------------------------------------------
  Attempts run:    ${#attempt_codes[@]} of $total_budget
  Total elapsed:   ${elapsed}s
  Final exit:      $final_exit ($( [[ $final_exit -eq 0 ]] && echo "SUCCESS" || echo "FAILED" ))
--------------------------------------------------------------------------------
  Per-attempt:
EOF
for i in "${!attempt_codes[@]}"; do
    n=$((i + 1))
    printf '    %d) exit=%-3s  %s\n' "$n" "${attempt_codes[$i]}" "${attempt_reasons[$i]}"
done
echo "================================================================================"

exit "$final_exit"
