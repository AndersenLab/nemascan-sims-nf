#!/bin/bash
#SBATCH --job-name=nemascan-resume
#SBATCH --account=eande106
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=nemascan-resume_%j.out
#SBATCH --error=nemascan-resume_%j.err
#
# Submit a resumed Nextflow run from the repository root on Rockfish.
#
# This job holds the Nextflow driver process. Child pipeline tasks are
# dispatched by Nextflow via its own SLURM executor; this job only runs
# the driver. --mem=64G lets the JVM auto-size its heap to ~16 GB via
# container cgroup support — sufficient for large runs (~200K tasks).
# Raise --mem for larger parameter grids.
#
# Usage (submit from the nemascan-sims-nf repo root on a login node):
#   sbatch scripts/resume_panelsize_ce.sbatch
#
# Do NOT submit from an interactive srun shell — Rockfish blocks sbatch
# from interactive srun allocations, which would prevent Nextflow's SLURM
# executor from dispatching child tasks.

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
NF_CONDA_ENV="${NF_CONDA_ENV:-/data/eande106/software/conda_envs/nf24_env}"
NXF_VER="${NXF_VER:-24.10.1}"

module load anaconda3/2022.05
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$NF_CONDA_ENV"
export NXF_VER

# Run from the directory where sbatch was invoked so .nextflow.log and the
# resume cache (.nextflow/) land in the repo root as expected.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "=== SLURM driver job ==="
echo "job_id:    ${SLURM_JOB_ID:-unknown}"
echo "node:      $(hostname)"
echo "cwd:       $PWD"
echo "NXF_VER:   $NXF_VER"
echo "conda env: $NF_CONDA_ENV"
echo

# ---------------------------------------------------------------------------
# Run parameters
# ---------------------------------------------------------------------------
VALIDATION_ROOT=/scratch4/eande106/Ryan/cme-test/2026-06-17

OUTDIR=${VALIDATION_ROOT}/cme_test_output
WORKDIR=${VALIDATION_ROOT}/cme_test_workdir

nextflow run main.nf \
    -resume \
    -profile rockfish \
    --strainfile data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/strains_panelsize_ce.tsv \
    --nqtl       data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/nqtl.csv \
    --h2         data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/h2.csv \
    --effect     data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/effect_sizes.csv \
    --reps       50 \
    --cv_maf     0.05 \
    --cv_ld      0.8 \
    --output_dir "${OUTDIR}" \
    -w           "${WORKDIR}" \
    -with-trace
