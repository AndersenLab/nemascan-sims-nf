#!/usr/bin/env bash
# compare_vp_reml.sh
#
# Run REML on the full GRM with the upscaled phenotype, then iteratively
# scale the phenotype by 1000x and re-run REML to find the scaling level
# at which GCTA produces a non-zero Vp estimate. Once found, run the
# fastGWA-mlm-exact mapping command at that scale.
#
# Usage:
#   cd issues/143-low-vp-error/data
#   bash compare_vp_reml.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/raw"
PROC_DIR="${SCRIPT_DIR}/proc"
CONTAINER="quay.io/biocontainers/gcta:1.94.1--h9ee0642_0"

MAKE_GRM_DIR="${RAW_DIR}/make_grm"
PERFORM_GWA_DIR="${RAW_DIR}/perform_gwa"
PHENO="${RAW_DIR}/check_vp/1_25_0.2_0.05_gamma_ce.full_sims.pheno"

BFILE_PREFIX="TO_SIMS_1_25_0.2_0.05_gamma_ce.full"

# Number of sequential 1000x scaling rounds (original + N scaled)
MAX_ROUNDS=4

# --- Preflight ----------------------------------------------------------------

if [ ! -f "$PHENO" ]; then
    echo "ERROR: $PHENO not found." >&2
    exit 1
fi

mkdir -p "$PROC_DIR"

# --- Run REML at each scaling level -------------------------------------------

CURRENT_PHENO="$PHENO"
CUMULATIVE_SCALE=1
NONZERO_LABEL=""
NONZERO_PHENO=""

for round in $(seq 0 "$MAX_ROUNDS"); do

    if [ "$round" -eq 0 ]; then
        LABEL="scale_1x"
        echo "=== REML: scale 1x (original upscaled phenotype) ==="
    else
        CUMULATIVE_SCALE=$((CUMULATIVE_SCALE * 1000))
        LABEL="scale_${CUMULATIVE_SCALE}x"
        echo ""
        echo "=== REML: scale ${CUMULATIVE_SCALE}x ==="

        # Scale phenotype values by 1000x
        SCALED_PHENO="${PROC_DIR}/pheno_${LABEL}.txt"
        awk '{printf "%s %s %.10g\n", $1, $2, $3 * 1000}' "$CURRENT_PHENO" > "$SCALED_PHENO"
        CURRENT_PHENO="$SCALED_PHENO"
    fi

    docker run --rm \
        -v "${MAKE_GRM_DIR}:/grm:ro" \
        -v "${CURRENT_PHENO}:/data/pheno.txt:ro" \
        -v "${PROC_DIR}:/out" \
        -w /work \
        "$CONTAINER" \
        bash -c '
            mkdir -p /work
            cp /grm/TO_SIMS_1_25_0.2_0.05_gamma_ce.full_gcta_grm_inbred.grm.bin   /work/full_grm.grm.bin
            cp /grm/TO_SIMS_1_25_0.2_0.05_gamma_ce.full_gcta_grm_inbred.grm.N.bin /work/full_grm.grm.N.bin
            cp /grm/TO_SIMS_1_25_0.2_0.05_gamma_ce.full_gcta_grm_inbred.grm.id    /work/full_grm.grm.id

            gcta64 --grm /work/full_grm \
                   --pheno /data/pheno.txt \
                   --reml \
                   --out /out/reml_'"$LABEL"' \
                   --thread-num 1

            echo "--- REML complete ---"
        '

    # Check if Vp is non-zero
    HSQ_FILE="${PROC_DIR}/reml_${LABEL}.hsq"
    VP=$(awk -F'\t' '$1 == "Vp" {print $2}' "$HSQ_FILE")
    if [ "$VP" != "0.000000" ] && [ -z "$NONZERO_LABEL" ]; then
        NONZERO_LABEL="$LABEL"
        NONZERO_PHENO="$CURRENT_PHENO"
        echo ""
        echo ">>> Vp transitioned to non-zero at ${LABEL} (Vp = ${VP})"
        break
    fi
done

# --- Summary ------------------------------------------------------------------

echo ""
echo "=============================================="
echo "  Vp COMPARISON ACROSS SCALING LEVELS"
echo "=============================================="
echo ""
printf "%-20s  %12s  %12s  %12s  %12s\n" "Scale" "V(G)" "V(e)" "Vp" "SE(Vp)"
printf "%-20s  %12s  %12s  %12s  %12s\n" "----" "----" "----" "----" "----"

for hsq_file in "$PROC_DIR"/reml_scale_*.hsq; do
    [ -e "$hsq_file" ] || continue
    label=$(basename "$hsq_file" .hsq | sed 's/reml_//')
    vp=$(awk -F'\t' '$1 == "Vp" {print $2}' "$hsq_file")
    vp_se=$(awk -F'\t' '$1 == "Vp" {print $3}' "$hsq_file")
    vg=$(awk -F'\t' '$1 == "V(G)" {print $2}' "$hsq_file")
    ve=$(awk -F'\t' '$1 == "V(e)" {print $2}' "$hsq_file")
    printf "%-20s  %12s  %12s  %12s  %12s\n" "$label" "$vg" "$ve" "$vp" "$vp_se"
done

# --- Run mapping at the non-zero Vp scale -------------------------------------

if [ -z "$NONZERO_LABEL" ]; then
    echo ""
    echo "WARNING: Vp remained zero across all scaling levels. No mapping run."
    exit 0
fi

echo ""
echo "=============================================="
echo "  MAPPING: fastGWA-mlm-exact at ${NONZERO_LABEL}"
echo "=============================================="
echo ""

# Check for required PLINK files
PLINK_DIR="${PERFORM_GWA_DIR}"
MISSING_FILES=()
for ext in bed bim fam; do
    if [ ! -f "${PLINK_DIR}/${BFILE_PREFIX}.${ext}" ]; then
        MISSING_FILES+=("${PLINK_DIR}/${BFILE_PREFIX}.${ext}")
    fi
done

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo "ERROR: PLINK files required for mapping not found:" >&2
    for f in "${MISSING_FILES[@]}"; do
        echo "  $f" >&2
    done
    echo "" >&2
    echo "Copy from HPC (GCTA_PERFORM_GWA workdir):" >&2
    echo "  rsync -avL rmckeow1@login.rockfish.jhu.edu:<workdir>/${BFILE_PREFIX}.{bed,bim,fam} ${PLINK_DIR}/" >&2
    exit 1
fi

docker run --rm \
    -v "${PLINK_DIR}:/plink:ro" \
    -v "${MAKE_GRM_DIR}:/full_grm:ro" \
    -v "${PERFORM_GWA_DIR}:/sparse:ro" \
    -v "${NONZERO_PHENO}:/data/pheno.txt:ro" \
    -v "${PROC_DIR}:/out" \
    -w /work \
    "$CONTAINER" \
    bash -c '
        mkdir -p /work

        # PLINK files
        cp /plink/'"$BFILE_PREFIX"'.bed /work/bfile.bed
        cp /plink/'"$BFILE_PREFIX"'.bim /work/bfile.bim
        cp /plink/'"$BFILE_PREFIX"'.fam /work/bfile.fam

        # Full GRM (needed by fastGWA-mlm-exact alongside sparse)
        cp /full_grm/TO_SIMS_1_25_0.2_0.05_gamma_ce.full_gcta_grm_inbred.grm.bin   /work/full_grm.grm.bin
        cp /full_grm/TO_SIMS_1_25_0.2_0.05_gamma_ce.full_gcta_grm_inbred.grm.N.bin /work/full_grm.grm.N.bin
        cp /full_grm/TO_SIMS_1_25_0.2_0.05_gamma_ce.full_gcta_grm_inbred.grm.id    /work/full_grm.grm.id

        # Sparse GRM
        cp /sparse/1_25_0.2_0.05_gamma_ce.full_sparse_grm_inbred.grm.sp /work/sparse_grm.grm.sp
        cp /sparse/1_25_0.2_0.05_gamma_ce.full_sparse_grm_inbred.grm.id /work/sparse_grm.grm.id

        # SNP list from .bim
        awk "{print \$2}" /work/bfile.bim > /work/plink_snplist.txt

        gcta64 --fastGWA-mlm-exact \
               --bfile /work/bfile \
               --grm-sparse /work/sparse_grm \
               --pheno /data/pheno.txt \
               --extract /work/plink_snplist.txt \
               --out /out/gwa_'"$NONZERO_LABEL"' \
               --thread-num 1

        echo "--- Mapping complete ---"
    '

echo ""
echo "Mapping output: ${PROC_DIR}/gwa_${NONZERO_LABEL}.fastGWA"
echo "Full .hsq files saved to: ${PROC_DIR}/"
