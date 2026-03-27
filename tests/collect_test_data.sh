#!/usr/bin/env bash
# Runs the test pipeline and collects outputs for integration tests.
#
# Usage:
#   bash tests/collect_test_data.sh [--profile <name>] [--no-legacy] [--collect-only] [--work-dir <path>] [--clean]
#
# Profiles:
#   test             - Fixed architecture: 1 population, 5 nQTL, 0.8 h², 1 rep (default)
#   test_variable    - Variable architecture: 8 nQTL × 4 h² × 2 reps = 256 mappings
#   test_three_species - 3 species × 3 populations = 9 groups, 36 mappings
#   test_cv_pool     - CV pool broader than marker set (cv_maf=0.01 < ms_maf=0.05)
#   test_hpc         - Fixed architecture, full Rockfish VCF (run on Rockfish with rockfish profile)
#   test_hpc_variable  - Variable architecture, full Rockfish VCF (run on Rockfish)
#   test_hpc_three_species - 3 species, full Rockfish VCFs (run on Rockfish)
#   test_hpc_cv_pool - CV pool, full Rockfish VCF (--no-legacy required; run on Rockfish)
#
# Options:
#   --no-legacy      Omit --legacy_assess (required for test_cv_pool: legacy path is
#                    incompatible when non-marker causal variants are present)
#   --collect-only   Skip pipeline run; collect from Analysis_Results-<profile>/ (or
#                    most recent Analysis_Results-* as fallback)
#   --work-dir <path> Collect cross-validation files (.fastGWA, .mlma, EIGEN) from this
#                    Nextflow work directory into tests/integration_data/<profile>/work/.
#                    For --collect-only with HPC profiles, pass the scratch work-dir path.
#   --clean          Remove previous profile artifacts and work directory before running
#
# Artifacts are stored in tests/integration_data/<profile>/ so each profile's output
# can be kept and tested independently.
#
# Requires: Docker for local runs, NXF_VER=24.10.4 (or any 24.10.x)
# Note: HPC profiles (test_hpc*) are run with -profile <profile>,rockfish directly on
# Rockfish, not via this script. Use --collect-only after the run to gather outputs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PROFILE="test"
NO_LEGACY=false
MODE="run"
DO_CLEAN=false
WORK_DIR_OVERRIDE=""

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
    --work-dir)
      WORK_DIR_OVERRIDE="$2"
      shift 2
      ;;
    --clean)
      DO_CLEAN=true
      shift
      ;;
    -h|--help)
      sed -n '2,23p' "$0"
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
  # HPC profiles reference VCFs on the Rockfish /vast filesystem — skip local check
  if [[ "$PROFILE" == test_hpc* ]]; then
    return 0
  fi
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
    local vcf="$PROJECT_DIR/data/test/test_ce.vcf.gz"
    if [[ ! -f "$vcf" ]]; then
      echo "ERROR: Missing VCF target: data/test/test_ce.vcf.gz" >&2
      echo "       Generate it with:" >&2
      echo "         bash data/test/generate_test_vcf.sh /path/to/WI.YYYYMMDD.hard-filter.isotype.vcf.gz" >&2
      echo "       See data/test/release_ids.sh for the canonical release date." >&2
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
    --output_dir "Analysis_Results-${PROFILE}" \
    -work-dir "$WORK_DIR"
fi

# --- Collect outputs ---
# Prefer the profile-named output directory (from --output_dir Analysis_Results-<profile>);
# fall back to the most recent date-stamped Analysis_Results-* for backward compatibility.
NAMED_DIR="$PROJECT_DIR/Analysis_Results-${PROFILE}"
if [[ -d "$NAMED_DIR" ]]; then
  RESULTS_DIR="$NAMED_DIR"
else
  RESULTS_DIR=$(ls -dt "$PROJECT_DIR"/Analysis_Results-* 2>/dev/null | head -1)
fi
if [[ -z "${RESULTS_DIR:-}" ]]; then
  echo "ERROR: No Analysis_Results-${PROFILE}/ or Analysis_Results-* directory found." >&2
  echo "       Run the pipeline with --output_dir Analysis_Results-${PROFILE} first." >&2
  exit 1
fi

echo "Collecting outputs from: $RESULTS_DIR"
mkdir -p "$OUTPUT_DIR"

# Copy DB directory (remove stale copy first to prevent nested db/db/)
if [[ -d "$RESULTS_DIR/db" ]]; then
  rm -rf "$OUTPUT_DIR/db"
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

# --- Collect cross-validation files from work directory ---
# Only the GWA output files (.fastGWA, .mlma) and EIGEN independent test files are needed
# for test-cross_validation.R. Collecting these selectively avoids syncing the full work
# directory (which can be 100+ GB on HPC) — the collected files are typically < 5% of that.
collect_cv_files() {
  local src_work="$1"
  local dest_work="$OUTPUT_DIR/work"

  if [[ ! -d "$src_work" ]]; then
    echo "WARNING: Work directory not found: $src_work — skipping cross-validation file collection" >&2
    return
  fi

  rm -rf "$dest_work"
  mkdir -p "$dest_work"

  local count=0

  # GWA output files (.fastGWA and .mlma) — flat copy, filenames contain all metadata
  while IFS= read -r -d '' f; do
    cp "$f" "$dest_work/"
    count=$((count + 1))
  done < <(find "$src_work" \( -name "*.fastGWA" -o -name "*.mlma" \) -print0)

  # EIGEN independent test count files
  while IFS= read -r -d '' f; do
    cp "$f" "$dest_work/"
    count=$((count + 1))
  done < <(find "$src_work" -name "*_total_independent_tests.txt" -print0)

  if [[ "$count" -eq 0 ]]; then
    echo "WARNING: No cross-validation files found in $src_work" >&2
    rm -rf "$dest_work"
  else
    echo "Collected $count cross-validation files to: $dest_work"
  fi
}

# Determine which work directory to collect from
CV_WORK_DIR="${WORK_DIR_OVERRIDE:-}"
if [[ -z "$CV_WORK_DIR" ]] && [[ "$MODE" != "collect" ]]; then
  # Local run: work directory is known
  CV_WORK_DIR="$WORK_DIR"
fi
if [[ -n "$CV_WORK_DIR" ]]; then
  collect_cv_files "$CV_WORK_DIR"
fi

echo ""
echo "Test data collected to: $OUTPUT_DIR"
echo ""
echo "Run integration tests with:"
echo ""

# Determine the expected population name for HPC vs local profiles
TEST_POP_FLAG=""
if [[ "$PROFILE" == test_hpc* ]]; then
  TEST_POP_FLAG="  TEST_POPULATION=ce.hpc \\"$'\n'
fi

# Include TEST_WORK_DIR if cross-validation files were collected
WORK_DIR_FLAG=""
if [[ -d "$OUTPUT_DIR/work" ]]; then
  WORK_DIR_FLAG="  TEST_WORK_DIR=%s/work \\\\\n"
fi

if [[ "$PROFILE" == *cv_pool* ]] || [[ "$NO_LEGACY" == "true" ]]; then
  # cv_pool: no legacy assessment; use TEST_CV_POOL flag for cv_pool-specific assertions
  CV_POOL_FLAG=""
  if [[ "$PROFILE" == *cv_pool* ]]; then
    CV_POOL_FLAG="  TEST_CV_POOL=true \\"$'\n'
  fi
  printf "  TEST_DB_DIR=%s/db \\\\\n" "$OUTPUT_DIR"
  printf "%s" "$TEST_POP_FLAG"
  if [[ -n "$WORK_DIR_FLAG" ]]; then printf "$WORK_DIR_FLAG" "$OUTPUT_DIR"; fi
  printf "  TEST_DB_ASSESSMENT=%s/db_simulation_assessment_results.tsv \\\\\n" "$OUTPUT_DIR"
  printf "%s" "$CV_POOL_FLAG"
  printf "  Rscript tests/run_tests.R\n"
else
  printf "  TEST_DB_DIR=%s/db \\\\\n" "$OUTPUT_DIR"
  printf "%s" "$TEST_POP_FLAG"
  if [[ -n "$WORK_DIR_FLAG" ]]; then printf "$WORK_DIR_FLAG" "$OUTPUT_DIR"; fi
  printf "  TEST_LEGACY_ASSESSMENT=%s/simulation_assessment_results.tsv \\\\\n" "$OUTPUT_DIR"
  printf "  TEST_DB_ASSESSMENT=%s/db_simulation_assessment_results.tsv \\\\\n" "$OUTPUT_DIR"
  printf "  Rscript tests/run_tests.R\n"
fi
