# bin/failure_trap.sh
#
# Sourced by every post-fanout process script after exporting the required
# env vars. Installs an EXIT trap that writes a JSON failure record into the
# shared .failures/ directory on every failed exit, and deletes that record
# on a successful exit. The file remaining in .failures/ at workflow end is
# therefore the most recent state — failed cells stay, retry-recovered cells
# are cleaned up.
#
# Why no maxRetries gate: dynamic errorStrategy callbacks (e.g. the rockfish
# config's `task.attempt <= N ? 'retry' : 'ignore'`) can return 'ignore'
# before maxRetries is reached. The previous gate `attempt > maxRetries`
# would silently skip those final-attempt failures. Per-cell keying with
# success-cleanup is independent of retry semantics.
#
# Required env vars (process script must export before sourcing):
#   NF_TRAP_FAILURES_DIR   — ${workflow.outputDir}/.failures
#   NF_TRAP_CELL_KEY       — process-built filename-safe string that uniquely
#                            identifies the cell at this DAG position. Must
#                            include every fanout variable in scope so that
#                            sibling tasks of this process never collide.
#   NF_TRAP_PAYLOAD        — process-built JSON body, MUST end in `}`. The
#                            trap appends an `,"exit":N` field by stripping
#                            the trailing `}` and re-closing.
#
# Forward-compat contract:
#   When a new fanout variable is added (e.g. species, future rep-disambiguators),
#   include it in BOTH NF_TRAP_CELL_KEY (so cell-key uniqueness scales) and
#   NF_TRAP_PAYLOAD (so replay.tsv gets it). NO edit to this script is needed.
#
# WARNING for test-injection:
#   Any synthetic failure injection (e.g. `exit 1` for Check 6) MUST be placed
#   AFTER the `source ... failure_trap.sh` line. Placing it earlier bypasses
#   trap registration and produces empty .failures/.

: "${NF_TRAP_FAILURES_DIR:?NF_TRAP_FAILURES_DIR must be set before sourcing failure_trap.sh}"
: "${NF_TRAP_CELL_KEY:?NF_TRAP_CELL_KEY must be set before sourcing failure_trap.sh}"
: "${NF_TRAP_PAYLOAD:?NF_TRAP_PAYLOAD must be set before sourcing failure_trap.sh}"

_handle_exit() {
    local rc=$?
    local final="${NF_TRAP_FAILURES_DIR}/${NF_TRAP_CELL_KEY}.json"

    if [ "$rc" -eq 0 ]; then
        rm -f "$final" 2>/dev/null || true
        return 0
    fi

    local tmp="${final}.tmp"
    # Splice exit code into payload: strip trailing `}` and append `,"exit":N}`.
    local body_sans_close="${NF_TRAP_PAYLOAD%\}}"
    printf '%s,"exit":%s}\n' "$body_sans_close" "$rc" \
        > "$tmp" 2>/dev/null \
        && mv "$tmp" "$final" 2>/dev/null \
        || true
}

trap _handle_exit EXIT
