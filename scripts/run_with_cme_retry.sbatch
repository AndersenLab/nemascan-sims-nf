#!/usr/bin/env bash
#SBATCH --job-name=nemascan-cme-retry
#SBATCH --account=eande106
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=nemascan-cme-retry_%j.out
#SBATCH --error=nemascan-cme-retry_%j.err
#
# SLURM wrapper for nemascan_run_with_cme_retry.sh on Rockfish.
#
# Submits a single SLURM job that runs the Nextflow driver under the
# auto-resume wrapper. The wrapper detects the issue #175 CME-at-706 abort
# and re-launches Nextflow with `-resume` up to --max-retries times. Child
# pipeline tasks are submitted by Nextflow via its own slurm executor as
# usual; this job only holds the driver.
#
# Memory (issue #175 OOM follow-up — issues/new-unexpected-error/rca.md):
#   The driver holds the whole task graph in its JVM heap for the entire run.
#   The original failure ran the driver on a login node where the per-user cgroup
#   capped RAM and the container-aware JVM sized its max heap to ~1.25 GB, so a
#   50-rep (~200K-task) run OOM'd. This job requests --mem=64G; the JVM then
#   auto-sizes its heap to ~25% of the cgroup (~16 GB) via container support — no
#   explicit -Xmx and no diagnostic flags are set, by design. Raise --mem for
#   larger runs.
#
# Submit from a login node with `sbatch`. Do NOT try to run the wrapper
# directly from an interactive `srun --pty` shell — Rockfish blocks
# sbatch submission from interactive srun allocations, which would
# prevent Nextflow's slurm executor from dispatching child tasks. The
# sbatch job below is a batch allocation, not an srun shell, so this
# restriction does not apply once the job starts.
#
# Usage:
#   sbatch scripts/run_with_cme_retry.sbatch \
#       -- main.nf -profile rockfish \
#          --strainfile strains.tsv \
#          --vcf 20231213 \
#          --output_dir /vast/eande106/${USER}/sim-$(date +%F)
#
# Override defaults (e.g. longer wall-clock, different retry budget):
#   sbatch --time=120:00:00 scripts/run_with_cme_retry.sbatch \
#       --max-retries 5 -- main.nf -profile rockfish ...
#
# The wrapper script accepts these flags BEFORE the `--`:
#   --max-retries N      (default 3)
#   --nextflow CMD       (default `nextflow` from PATH)
#   --log-file PATH      (default ./.nextflow.log)
#   --no-resume-first    (default: -resume IS passed on attempt 1)
#
# Everything AFTER `--` is forwarded verbatim to `nextflow run`.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    cat >&2 <<USAGE
usage: sbatch $0 [wrapper-opts] -- <args to 'nextflow run'>

  wrapper-opts:
    --max-retries N        (default 3)
    --nextflow CMD         (default 'nextflow')
    --log-file PATH        (default ./.nextflow.log)
    --no-resume-first

  example:
    sbatch $0 --max-retries 3 -- main.nf -profile rockfish \\
        --strainfile strains.tsv --vcf 20231213 \\
        --output_dir /vast/eande106/\${USER}/sim-\$(date +%F)
USAGE
    exit 1
fi

# Resolve repo root from the script location so the user can submit from
# any cwd. SLURM sets $SLURM_SUBMIT_DIR to the caller's cwd; we honor that
# so Nextflow writes .nextflow.log and resumes from the same place across
# retries.
wrapper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$wrapper_dir/.." && pwd)"
wrapper="$wrapper_dir/nemascan_run_with_cme_retry.sh"

if [[ ! -x "$wrapper" ]]; then
    echo "ERR: wrapper not executable at $wrapper" >&2
    exit 1
fi

# Rockfish nextflow conda env — same one the pipeline docs and CLAUDE.md
# call out.
NF_CONDA_ENV="${NF_CONDA_ENV:-/data/eande106/software/conda_envs/nf24_env}"
NXF_VER="${NXF_VER:-24.10.1}"

# Activate conda. `conda activate` requires the shell hook to be sourced;
# `source activate` is the portable fallback on managed envs.
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$NF_CONDA_ENV"
else
    # shellcheck disable=SC1091
    source activate "$NF_CONDA_ENV"
fi
export NXF_VER

# Run from $SLURM_SUBMIT_DIR (where the user invoked sbatch) so .nextflow.log
# and the resume cache land where the user expects.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "=== SLURM driver job ==="
echo "job_id:        ${SLURM_JOB_ID:-unknown}"
echo "node:          $(hostname)"
echo "submit_dir:    ${SLURM_SUBMIT_DIR:-$PWD}"
echo "repo_root:     $repo_root"
echo "wrapper:       $wrapper"
echo "NF_CONDA_ENV:  $NF_CONDA_ENV"
echo "NXF_VER:       $NXF_VER"
echo "args:          $*"
echo

exec "$wrapper" "$@"
