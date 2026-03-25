#!/usr/bin/env bash
#
# Generate the C. briggsae test VCF for nemascan-sims-nf integration testing.
#
# Subsets a CaeNDR C. briggsae isotype VCF to a small set of strains and
# chromosomes I, II, V. Removes monomorphic sites and strips contig headers
# for absent chromosomes (III, IV, X). Creates a dated symlink so that
# extractVcfReleaseId() in main.nf can resolve the 8-digit release date.
#
# Prerequisites:
#   - bcftools (>= 1.16)
#   - tabix (from htslib)
#   - Source VCF: CB.YYYYMMDD.hard-filter.isotype.vcf.gz
#     Download from CaeNDR: https://caendr.org/data/release/latest
#
# Usage:
#   ./generate_test_cb_vcf.sh <source_cb_vcf> [strain_list_file]
#
#   source_cb_vcf   Path to CB.YYYYMMDD.hard-filter.isotype.vcf.gz
#   strain_list_file  Optional: one strain name per line.
#                     Default: first 20 strains listed in the source VCF.
#
# Output:
#   data/test/test_cb.vcf.gz             (BGZF-compressed VCF)
#   data/test/test_cb.vcf.gz.tbi         (tabix index)
#   data/test/test_cb.YYYYMMDD.vcf.gz    (dated symlink → test_cb.vcf.gz)
#   data/test/test_cb.YYYYMMDD.vcf.gz.tbi (dated symlink → test_cb.vcf.gz.tbi)
#
# On completion the script writes the three C. briggsae rows directly into
# data/test/test_strains_three_species.txt, replacing any existing cb.* rows.
# Rows for other species (ce.*, ct.*) are preserved unchanged.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_VCF="${SCRIPT_DIR}/test_cb.vcf.gz"
CHROMOSOMES="I,II,V"
DEFAULT_STRAIN_COUNT=200
# POP_SPLIT_COUNT is calculated dynamically from actual strain count (see below)

# Fraction of variants to randomly retain. Values < 1.0 reduce PLINK marker counts,
# keeping EIGEN matrices small enough for local Docker execution.
# Set VARIANT_SAMPLE_FRACTION=1.0 to disable sampling (production use).
VARIANT_SAMPLE_FRACTION="${VARIANT_SAMPLE_FRACTION:-0.15}"
VARIANT_SAMPLE_SEED="${VARIANT_SAMPLE_SEED:-42}"  # RNG seed for reproducibility

# --- Parse arguments ---
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <source_cb_vcf> [strain_list_file]"
    echo ""
    echo "  source_cb_vcf   Path to CB.YYYYMMDD.hard-filter.isotype.vcf.gz"
    echo "  strain_list_file  Optional: one strain per line (default: first 20 in VCF)"
    exit 1
fi

SOURCE_VCF="$1"
STRAIN_LIST_FILE_ARG="${2:-}"

if [[ "$SOURCE_VCF" != http* ]] && [[ ! -f "$SOURCE_VCF" ]]; then
    echo "Error: Source VCF not found: $SOURCE_VCF"
    exit 1
fi

# --- Check dependencies ---
for cmd in bcftools tabix; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: $cmd not found on PATH"
        exit 1
    fi
done

# --- Extract release date from source VCF filename ---
VCF_BASENAME="$(basename "$SOURCE_VCF")"
DATE=$(echo "$VCF_BASENAME" | grep -oE '[0-9]{8}' | head -1)
if [[ -z "$DATE" ]]; then
    echo "Error: Cannot extract 8-digit release date from filename: $VCF_BASENAME"
    echo "Rename the file to include a YYYYMMDD date (e.g. CB.20221202.hard-filter.isotype.vcf.gz)"
    exit 1
fi
echo "Detected release date: $DATE"

# --- Build strain list ---
STRAIN_LIST_FILE=$(mktemp)
if [[ -n "$STRAIN_LIST_FILE_ARG" ]]; then
    if [[ ! -f "$STRAIN_LIST_FILE_ARG" ]]; then
        echo "Error: Strain list file not found: $STRAIN_LIST_FILE_ARG"
        exit 1
    fi
    cp "$STRAIN_LIST_FILE_ARG" "$STRAIN_LIST_FILE"
    echo "Using provided strain list: $STRAIN_LIST_FILE_ARG"
else
    echo "No strain list provided — using first ${DEFAULT_STRAIN_COUNT} strains from VCF header..."
    bcftools query -l "$SOURCE_VCF" | head -${DEFAULT_STRAIN_COUNT} > "$STRAIN_LIST_FILE"
fi

STRAIN_COUNT=$(wc -l < "$STRAIN_LIST_FILE" | tr -d ' ')
POP_SPLIT_COUNT=$(( (STRAIN_COUNT + 1) / 2 ))
echo "Subsetting ${STRAIN_COUNT} strains (popA=${POP_SPLIT_COUNT}, popB=$(( STRAIN_COUNT - POP_SPLIT_COUNT ))) to chromosomes ${CHROMOSOMES}..."

# --- Subset VCF ---
echo "Running bcftools view (this may take several minutes for a large source VCF)..."
TEMP_VCF=$(mktemp "${SCRIPT_DIR}/test_cb_XXXXXX.vcf.gz")
bcftools view \
    -S "$STRAIN_LIST_FILE" \
    -t "$CHROMOSOMES" \
    --min-ac 1 \
    -Oz -o "$TEMP_VCF" \
    "$SOURCE_VCF"

# --- Random variant sampling ---
# Randomly retain VARIANT_SAMPLE_FRACTION of variants to keep PLINK marker counts
# manageable for local Docker execution (target: ~15K–50K variants/chrom → ~2K–5K markers).
echo "Random variant sampling: keeping ${VARIANT_SAMPLE_FRACTION} fraction (seed ${VARIANT_SAMPLE_SEED})..."
SAMPLED_VCF=$(mktemp "${SCRIPT_DIR}/test_cb_sampled_XXXXXX.vcf.gz")
bcftools view "$TEMP_VCF" | \
    awk -v rate="${VARIANT_SAMPLE_FRACTION}" -v seed="${VARIANT_SAMPLE_SEED}" \
        'BEGIN{srand(seed)} /^#/{print; next} rand() < rate' | \
    bcftools view -Oz -o "$SAMPLED_VCF" -
rm -f "$TEMP_VCF"
TEMP_VCF="$SAMPLED_VCF"

# --- Clean contig headers ---
# Strip ##contig lines for chromosomes not in the subset (III, IV, X)
echo "Cleaning contig headers..."
HEADER_FILE=$(mktemp)
bcftools view -h "$TEMP_VCF" \
    | grep -v '##contig=<ID=III' \
    | grep -v '##contig=<ID=IV' \
    | grep -v '##contig=<ID=X,' \
    > "$HEADER_FILE"

bcftools reheader -h "$HEADER_FILE" -o "$OUTPUT_VCF" "$TEMP_VCF"

# --- Index ---
echo "Creating tabix index..."
tabix -p vcf "$OUTPUT_VCF"

# --- Create dated symlinks ---
SYMLINK_VCF="${SCRIPT_DIR}/test_cb.${DATE}.vcf.gz"
SYMLINK_TBI="${SCRIPT_DIR}/test_cb.${DATE}.vcf.gz.tbi"

# Remove existing symlinks with this date if they exist
rm -f "$SYMLINK_VCF" "$SYMLINK_TBI"

ln -s test_cb.vcf.gz     "$SYMLINK_VCF"
ln -s test_cb.vcf.gz.tbi "$SYMLINK_TBI"
echo "Created symlinks:"
echo "  $SYMLINK_VCF → test_cb.vcf.gz"
echo "  $SYMLINK_TBI → test_cb.vcf.gz.tbi"

# --- Cleanup ---
rm -f "$TEMP_VCF" "$STRAIN_LIST_FILE" "$HEADER_FILE"

# --- Verify ---
SAMPLE_COUNT=$(bcftools query -l "$OUTPUT_VCF" | wc -l | tr -d ' ')
VARIANT_COUNT=$(bcftools index -n "$OUTPUT_VCF")
FILE_SIZE=$(ls -lh "$OUTPUT_VCF" | awk '{print $5}')
CONTIGS=$(bcftools view -h "$OUTPUT_VCF" | grep '##contig' | sed 's/.*ID=//;s/,.*//')

echo ""
echo "=== C. briggsae test VCF generated ==="
echo "  File:       $OUTPUT_VCF"
echo "  Size:       $FILE_SIZE"
echo "  Samples:    $SAMPLE_COUNT"
echo "  Variants:   $VARIANT_COUNT"
echo "  Contigs:    $(echo $CONTIGS | tr '\n' ' ')"
echo ""

# --- Update test_strains_three_species.txt ---
THREE_SPECIES_FILE="${SCRIPT_DIR}/test_strains_three_species.txt"
FINAL_STRAIN_COUNT=$(bcftools query -l "$OUTPUT_VCF" | wc -l | tr -d ' ')
POPB_SPLIT_COUNT=$(( FINAL_STRAIN_COUNT - POP_SPLIT_COUNT ))
STRAINS_CSV=$(bcftools query -l "$OUTPUT_VCF" | tr '\n' ',' | sed 's/,$//')
POPA_CSV=$(bcftools query -l "$OUTPUT_VCF" | head -${POP_SPLIT_COUNT} | tr '\n' ',' | sed 's/,$//')
POPB_CSV=$(bcftools query -l "$OUTPUT_VCF" | tail -${POPB_SPLIT_COUNT} | tr '\n' ',' | sed 's/,$//')

{
    printf "group\tspecies\tvcf\tms_maf\tms_ld\tstrains\n"
    # Preserve non-cb rows from existing strainfile
    if [[ -f "$THREE_SPECIES_FILE" ]]; then
        awk -F'\t' 'NR>1 && $1 !~ /^cb\./' "$THREE_SPECIES_FILE"
    fi
    printf "cb.test\tc_briggsae\tdata/test/test_cb.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$DATE" "$STRAINS_CSV"
    printf "cb.test.popA\tc_briggsae\tdata/test/test_cb.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$DATE" "$POPA_CSV"
    printf "cb.test.popB\tc_briggsae\tdata/test/test_cb.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$DATE" "$POPB_CSV"
} > "${THREE_SPECIES_FILE}.tmp"
mv "${THREE_SPECIES_FILE}.tmp" "$THREE_SPECIES_FILE"
echo "Updated $THREE_SPECIES_FILE (cb.test, cb.test.popA, cb.test.popB)"

echo "Remember to:"
echo "  1. git add data/test/test_cb.${DATE}.vcf.gz data/test/test_cb.${DATE}.vcf.gz.tbi"
echo "  2. git add data/test/test_strains_three_species.txt"
echo "  3. Update CB_VCF_RELEASE in data/test/release_ids.sh to ${DATE}"
echo ""
echo "Done."
