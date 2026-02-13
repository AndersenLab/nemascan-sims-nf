#!/usr/bin/env bash
# Runs the test pipeline and collects outputs for integration tests.
#
# Usage:
#   bash tests/collect_test_data.sh              # Run pipeline + collect outputs
#   bash tests/collect_test_data.sh --collect-only  # Collect from previous run
#   bash tests/collect_test_data.sh --clean         # Remove previous results first
#
# Requires: Docker, NXF_VER=24.10.4 (or any 24.10.x)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$PROJECT_DIR/tests/integration_data"
WORK_DIR="$PROJECT_DIR/tests/.nf-work"

MODE="run"
for arg in "$@"; do
  case "$arg" in
    --collect-only) MODE="collect" ;;
    --clean)
      echo "Cleaning previous results..."
      rm -rf "$OUTPUT_DIR" "$WORK_DIR"
      # Remove any Analysis_Results directories
      rm -rf "$PROJECT_DIR"/Analysis_Results-*
      ;;
    -h|--help)
      sed -n '2,9p' "$0"
      exit 0
      ;;
    *) echo "Unknown option: $arg"; exit 1 ;;
  esac
done

# --- Run pipeline (unless --collect-only) ---
if [[ "$MODE" != "collect" ]]; then
  echo "Running test pipeline with --legacy_assess (both DB and legacy paths)..."
  NXF_VER="${NXF_VER:-24.10.4}" nextflow run "$PROJECT_DIR/main.nf" \
    -profile test,docker --legacy_assess \
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
echo "  TEST_DB_DIR=$OUTPUT_DIR/db \\"
echo "  TEST_WORK_DIR=$WORK_DIR \\"
echo "  TEST_EXISTING_ASSESSMENT=$OUTPUT_DIR/simulation_assessment_results.tsv \\"
echo "  TEST_DB_ASSESSMENT=$OUTPUT_DIR/db_simulation_assessment_results.tsv \\"
echo "  Rscript tests/run_tests.R"
