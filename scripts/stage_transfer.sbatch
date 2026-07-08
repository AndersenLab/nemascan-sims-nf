#!/usr/bin/env bash
#SBATCH --job-name=stage-transfer
#SBATCH --account=eande106
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=stage_transfer_%j.out
#SBATCH --error=stage_transfer_%j.err
#
# SLURM wrapper for pretransfer_manifest.sh.
#
# The staging step (tar --sort=name | zstd -T0 -9, plus per-file SHA-256 of
# an ~11 GB output tree) pins every allocated core for several minutes and
# does heavy I/O against /vast. That profile is not appropriate for a
# Rockfish login node — run it on a compute node instead via this script.
#
# Usage:
#   sbatch scripts/stage_transfer.sbatch <output_dir> [staging_dir] [--force]
#
# Example:
#   sbatch scripts/stage_transfer.sbatch \
#     /vast/eande106/$USER/Sims_ce-cb-ct_nqtl1-50_h2grid_50reps
#
# Override defaults on the command line if needed:
#   sbatch --account=myaccount --cpus-per-task=16 --time=2:00:00 \
#     scripts/stage_transfer.sbatch /vast/eande106/$USER/Sims_...
#
# Output/error logs land in stage_transfer_<jobid>.out / .err in the
# submission directory. The script at the end prints the same rsync +
# sha256sum -c + tar extract commands the user should run from their
# local machine.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: sbatch $0 <output_dir> [staging_dir] [--force]" >&2
    exit 1
fi

# Under sbatch, the script runs from a node-local spool copy, so BASH_SOURCE
# does not point at the repo. Use SLURM_SUBMIT_DIR (the submission cwd) when set.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    repo_root="$SLURM_SUBMIT_DIR"
else
    wrapper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    repo_root="$(cd "$wrapper_dir/.." && pwd)"
fi
manifest_tool="$repo_root/pretransfer_manifest.sh"

if [[ ! -x "$manifest_tool" ]]; then
    echo "ERR: pretransfer_manifest.sh not found or not executable at $manifest_tool" >&2
    exit 1
fi

echo "=== SLURM staging job ==="
echo "job_id:        ${SLURM_JOB_ID:-unknown}"
echo "node:          $(hostname)"
echo "cpus:          ${SLURM_CPUS_PER_TASK:-?}"
echo "submit_dir:    ${SLURM_SUBMIT_DIR:-$PWD}"
echo "repo_root:     $repo_root"
echo "manifest_tool: $manifest_tool"
echo "args:          $*"
echo

exec "$manifest_tool" "$@"
