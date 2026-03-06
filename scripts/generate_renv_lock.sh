#!/usr/bin/env bash
# generate_renv_lock.sh — regenerate renv.lock from the r_packages container
#
# Usage:
#   bash scripts/generate_renv_lock.sh docker      # local execution with Docker
#   bash scripts/generate_renv_lock.sh singularity # HPC execution with Singularity
#
# Writes renv.lock to the project root. Commit the result.

set -euo pipefail

CONTAINER="andersenlab/r_packages:20250519"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# --- Require explicit runtime argument ---
if [ $# -ne 1 ] || { [ "$1" != "docker" ] && [ "$1" != "singularity" ]; }; then
  echo "Usage: bash scripts/generate_renv_lock.sh <docker|singularity>"
  echo ""
  echo "  docker      — local execution via Docker (requires docker)"
  echo "  singularity — HPC execution via Singularity (Rockfish)"
  exit 1
fi

RUNTIME="$1"

echo "==> Generating renv.lock from ${CONTAINER} (${RUNTIME})"
echo "==> Project root: ${PROJECT_ROOT}"

# Shared renv snapshot expression
RSCRIPT_EXPR="
  renv::snapshot(
    project  = '/project',
    type     = 'all',
    lockfile = '/project/renv.lock',
    prompt   = FALSE
  )
  cat('renv.lock written to project root.\n')
"

if [ "${RUNTIME}" = "docker" ]; then
  if ! command -v docker &>/dev/null; then
    echo "ERROR: docker not found on PATH."
    exit 1
  fi
  echo "==> Pulling ${CONTAINER} ..."
  docker pull "${CONTAINER}"

  # Verify renv is installed in the container
  RENV_PRESENT=$(docker run --rm "${CONTAINER}" \
    Rscript -e "cat(nchar(system.file(package='renv')) > 0)" 2>/dev/null || echo "FALSE")
  if [ "${RENV_PRESENT}" != "TRUE" ]; then
    echo "ERROR: renv is not installed in ${CONTAINER}."
    echo "Install renv into the container image and rebuild before running this script."
    exit 1
  fi

  docker run --rm \
    --volume "${PROJECT_ROOT}:/project" \
    "${CONTAINER}" \
    Rscript -e "${RSCRIPT_EXPR}"

elif [ "${RUNTIME}" = "singularity" ]; then
  if ! command -v singularity &>/dev/null; then
    echo "ERROR: singularity not found on PATH."
    exit 1
  fi
  singularity exec \
    --bind "${PROJECT_ROOT}:/project" \
    "docker://${CONTAINER}" \
    Rscript -e "${RSCRIPT_EXPR}"
fi

echo "==> Done. Commit renv.lock to the repository."
