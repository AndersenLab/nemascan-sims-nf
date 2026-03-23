#!/usr/bin/env bash
# Runs the test pipeline and collects outputs for integration tests.
#
# Usage:
#   bash tests/collect_test_data.sh [--profile <name>] [--no-legacy] [--collect-only] [--clean]
#
# Profiles:
#   test             - Fixed architecture: 1 population, 5 nQTL, 0.8 h², 1 rep (default)
#   test_variable    - Variable architecture: 8 nQTL × 4 h² × 2 reps = 256 mappings
#   test_three_species - 3 species × 3 populations = 9 groups, 36 mappings
#   test_cv_pool     - CV pool broader than marker set (cv_maf=0.01 < ms_maf=0.05)
#
# Options:
#   --no-legacy      Omit --legacy_assess (required for test_cv_pool: legacy path is
#                    incompatible when non-marker causal variants are present)
#   --collect-only   Skip pipeline run; collect from most recent Analysis_Results-*
#   --clean          Remove previous profile artifacts and work directory before running
#
# Artifacts are stored in tests/integration_data/<profile>/ so each profile's output
# can be kept and tested independently.
#
# Requires: Docker, NXF_VER=24.10.4 (or any 24.10.x)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PROFILE="test"
NO_LEGACY=false
MODE="run"
DO_CLEAN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --no-legacy)
      NO_LEGACY=true
      shift
      ;;
    --collect-only)
      MODE="collect"
      shift
      ;;
    --clean)
      DO_CLEAN=true
      shift
      ;;
    -h|--help)
      sed -n '2,21p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

OUTPUT_DIR="$PROJECT_DIR/tests/integration_data/$PROFILE"
WORK_DIR="$PROJECT_DIR/tests/.nf-work"

if [[ "$DO_CLEAN" == "true" ]]; then
  echo "Cleaning previous results for profile: $PROFILE..."
  rm -rf "$OUTPUT_DIR" "$WORK_DIR"
  rm -rf "$PROJECT_DIR"/Analysis_Results-*
fi

# --- Pre-flight: verify VCF symlink targets exist ---
check_vcfs() {
  local missing=false
  if [[ "$PROFILE" == "test_three_species" ]]; then
    for spec in ce cb ct; do
      local vcf="$PROJECT_DIR/data/test/test_${spec}.vcf.gz"
      if [[ ! -f "$vcf" ]]; then
        echo "ERROR: Missing VCF target: data/test/test_${spec}.vcf.gz" >&2
        echo "       Generate it with:" >&2
        echo "         bash data/test/generate_test_${spec}_vcf.sh /path/to/source.vcf.gz" >&2
        missing=true
      fi
    done
  else
    local vcf="$PROJECT_DIR/data/test/test.vcf.gz"
    if [[ ! -f "$vcf" ]]; then
      echo "ERROR: Missing VCF target: data/test/test.vcf.gz" >&2
      echo "       Generate it with:" >&2
      echo "         bash data/test/generate_test_vcf.sh /path/to/WI.20220216.hard-filter.isotype.vcf.gz" >&2
      missing=true
    fi
  fi
  if [[ "$missing" == "true" ]]; then
    echo "" >&2
    echo "VCF files are gitignored and must be generated locally before running collect_test_data.sh." >&2
    exit 1
  fi
}

# --- Run pipeline (unless --collect-only) ---
if [[ "$MODE" != "collect" ]]; then
  check_vcfs
  LEGACY_FLAG=""
  if [[ "$NO_LEGACY" == "false" ]]; then
    LEGACY_FLAG="--legacy_assess"
  fi

  echo "Running test pipeline (profile: ${PROFILE}${NO_LEGACY:+, no legacy assess})..."
  NXF_VER="${NXF_VER:-24.10.4}" nextflow run "$PROJECT_DIR/main.nf" \
    -profile "${PROFILE},docker" $LEGACY_FLAG \
    -work-dir "$WORK_DIR"
fi

# --- Collect outputs ---
RESULTS_DIR=$(ls -dt "$PROJECT_DIR"/Analysis_Results-* 2>/dev/null | head -1)
if [[ -z "$RESULTS_DIR" ]]; then
  echo "ERROR: No Analysis_Results-* directory found. Run the pipeline first." >&2
  exit 1
fi

echo "Collecting outputs from: $RESULTS_DIR"
mkdir -p "$OUTPUT_DIR"

# Copy DB directory
if [[ -d "$RESULTS_DIR/db" ]]; then
  cp -R "$RESULTS_DIR/db" "$OUTPUT_DIR/db"
else
  echo "WARNING: No db/ directory found in $RESULTS_DIR" >&2
fi

# Copy assessment TSVs
for tsv in db_simulation_assessment_results.tsv simulation_assessment_results.tsv; do
  if [[ -f "$RESULTS_DIR/$tsv" ]]; then
    cp "$RESULTS_DIR/$tsv" "$OUTPUT_DIR/"
  else
    echo "WARNING: $tsv not found in $RESULTS_DIR" >&2
  fi
done

echo ""
echo "Test data collected to: $OUTPUT_DIR"
echo ""
echo "Run integration tests with:"
echo ""

if [[ "$PROFILE" == "test_cv_pool" ]] || [[ "$NO_LEGACY" == "true" ]]; then
  # cv_pool: no legacy assessment; use TEST_CV_POOL flag for cv_pool-specific assertions
  CV_POOL_FLAG=""
  if [[ "$PROFILE" == "test_cv_pool" ]]; then
    CV_POOL_FLAG="  TEST_CV_POOL=true \\"$'\n'
  fi
  printf "  TEST_DB_DIR=%s/db \\\\\n" "$OUTPUT_DIR"
  printf "  TEST_WORK_DIR=%s \\\\\n" "$WORK_DIR"
  printf "  TEST_DB_ASSESSMENT=%s/db_simulation_assessment_results.tsv \\\\\n" "$OUTPUT_DIR"
  printf "%s" "$CV_POOL_FLAG"
  printf "  Rscript tests/run_tests.R\n"
else
  printf "  TEST_DB_DIR=%s/db \\\\\n" "$OUTPUT_DIR"
  printf "  TEST_WORK_DIR=%s \\\\\n" "$WORK_DIR"
  printf "  TEST_LEGACY_ASSESSMENT=%s/simulation_assessment_results.tsv \\\\\n" "$OUTPUT_DIR"
  printf "  TEST_DB_ASSESSMENT=%s/db_simulation_assessment_results.tsv \\\\\n" "$OUTPUT_DIR"
  printf "  Rscript tests/run_tests.R\n"
fi
