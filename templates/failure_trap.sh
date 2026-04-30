# templates/failure_trap.sh
#
# Sourced by every post-fanout process script after exporting the required
# tuple env vars. Installs an EXIT trap that, on the FINAL retry attempt only,
# atomically writes a JSON failure record into the shared .failures/ directory.
#
# Required env vars (process script must export before sourcing):
#   NF_TRAP_SESSION_ID    — workflow.sessionId, baked at script-render time
#   NF_TRAP_FAILURES_DIR  — ${workflow.outputDir}/.failures
#   NF_TRAP_TASK_HASH     — ${task.hash}
#   NF_TRAP_ATTEMPT       — ${task.attempt}
#   NF_TRAP_MAX_RETRIES   — ${task.maxRetries}
#   GROUP, MAF, NQTL, EFFECT, H2, REP, MODE, TYPE — eight-tuple

: "${NF_TRAP_SESSION_ID:?NF_TRAP_SESSION_ID must be set before sourcing failure_trap.sh}"
: "${NF_TRAP_FAILURES_DIR:?NF_TRAP_FAILURES_DIR must be set before sourcing failure_trap.sh}"
: "${NF_TRAP_TASK_HASH:?NF_TRAP_TASK_HASH must be set before sourcing failure_trap.sh}"
: "${NF_TRAP_ATTEMPT:?NF_TRAP_ATTEMPT must be set before sourcing failure_trap.sh}"
: "${NF_TRAP_MAX_RETRIES:?NF_TRAP_MAX_RETRIES must be set before sourcing failure_trap.sh}"
: "${GROUP:?GROUP must be set before sourcing failure_trap.sh}"
: "${MAF:?MAF must be set before sourcing failure_trap.sh}"
: "${NQTL:?NQTL must be set before sourcing failure_trap.sh}"
: "${EFFECT:?EFFECT must be set before sourcing failure_trap.sh}"
: "${H2:?H2 must be set before sourcing failure_trap.sh}"
: "${REP:?REP must be set before sourcing failure_trap.sh}"
: "${MODE:?MODE must be set before sourcing failure_trap.sh}"
: "${TYPE:?TYPE must be set before sourcing failure_trap.sh}"

_record_failure() {
    local rc=$?
    [ "$rc" -ne 0 ] || return 0
    [ "$NF_TRAP_ATTEMPT" -gt "$NF_TRAP_MAX_RETRIES" ] || return 0

    local tmp="${NF_TRAP_FAILURES_DIR}/${NF_TRAP_TASK_HASH}.json.tmp"
    local final="${NF_TRAP_FAILURES_DIR}/${NF_TRAP_TASK_HASH}.json"

    printf '{"session":"%s","group":"%s","maf":%s,"nqtl":"%s","effect":"%s","h2":%s,"rep":%s,"mode":"%s","type":"%s","attempt":%s,"exit":%s}\n' \
        "$NF_TRAP_SESSION_ID" \
        "$GROUP" "$MAF" "$NQTL" "$EFFECT" "$H2" "$REP" "$MODE" "$TYPE" \
        "$NF_TRAP_ATTEMPT" "$rc" \
        > "$tmp" 2>/dev/null \
        && mv "$tmp" "$final" 2>/dev/null \
        || true
}

trap _record_failure EXIT
