#!/usr/bin/env bash
#
# Generate the C. elegans test VCF for nemascan-sims-nf integration testing.
#
# Subsets a CaeNDR C. elegans isotype VCF to a set of strains and
# chromosomes I, II, V. Removes monomorphic sites and strips contig headers
# for absent chromosomes (III, IV, X). Creates a dated symlink so that
# extractVcfReleaseId() in main.nf can resolve the 8-digit release date.
#
# Prerequisites:
#   - bcftools (>= 1.16)
#   - tabix (from htslib)
#   - Source VCF: WI.YYYYMMDD.hard-filter.isotype.vcf.gz
#     Download from CaeNDR:
#       https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/c_elegans/YYYYMMDD/variation/WI.YYYYMMDD.hard-filter.isotype.vcf.gz
#
# Usage:
#   ./generate_test_ce_vcf.sh <source_ce_vcf> [strain_list_file]
#
#   source_ce_vcf     Path to WI.YYYYMMDD.hard-filter.isotype.vcf.gz or an HTTPS URL
#   strain_list_file  Optional: one strain name per line.
#                     Default: first 200 strains listed in the source VCF.
#
# Output:
#   data/test/test_ce.vcf.gz             (BGZF-compressed VCF)
#   data/test/test_ce.vcf.gz.tbi         (tabix index)
#   data/test/test_ce.YYYYMMDD.vcf.gz    (dated symlink → test_ce.vcf.gz)
#   data/test/test_ce.YYYYMMDD.vcf.gz.tbi (dated symlink → test_ce.vcf.gz.tbi)
#
# On completion the script writes the three C. elegans rows directly into
# data/test/test_strains_three_species.txt, replacing any existing ce.* rows.
# Rows for other species (cb.*, ct.*) are preserved unchanged.
#
# Note: This script outputs to test_ce.vcf.gz. generate_test_vcf.sh is a
# convenience wrapper that reads strains from the source VCF header and also
# outputs to test_ce.vcf.gz with a dated symlink, and rewrites test_strains.txt.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_VCF="${SCRIPT_DIR}/test_ce.vcf.gz"
CHROMOSOMES="I,II,V"
DEFAULT_STRAIN_COUNT=200
POP_SPLIT_COUNT=100   # strains per popA / popB (must be <= DEFAULT_STRAIN_COUNT / 2)

# --- Parse arguments ---
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <source_ce_vcf> [strain_list_file]"
    echo ""
    echo "  source_ce_vcf     Path to WI.YYYYMMDD.hard-filter.isotype.vcf.gz or HTTPS URL"
    echo "  strain_list_file  Optional: one strain per line (default: first ${DEFAULT_STRAIN_COUNT} in VCF)"
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
    echo "Rename the file to include a YYYYMMDD date (e.g. WI.20250625.hard-filter.isotype.vcf.gz)"
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
echo "Subsetting ${STRAIN_COUNT} strains to chromosomes ${CHROMOSOMES}..."

# --- Subset VCF ---
echo "Running bcftools view (this may take several minutes for a large source VCF)..."
TEMP_VCF=$(mktemp "${SCRIPT_DIR}/test_ce_XXXXXX.vcf.gz")
bcftools view \
    -S "$STRAIN_LIST_FILE" \
    -t "$CHROMOSOMES" \
    --min-ac 1 \
    -Oz -o "$TEMP_VCF" \
    "$SOURCE_VCF"

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
SYMLINK_VCF="${SCRIPT_DIR}/test_ce.${DATE}.vcf.gz"
SYMLINK_TBI="${SCRIPT_DIR}/test_ce.${DATE}.vcf.gz.tbi"

# Remove existing symlinks with this date if they exist
rm -f "$SYMLINK_VCF" "$SYMLINK_TBI"

ln -s test_ce.vcf.gz     "$SYMLINK_VCF"
ln -s test_ce.vcf.gz.tbi "$SYMLINK_TBI"
echo "Created symlinks:"
echo "  $SYMLINK_VCF → test_ce.vcf.gz"
echo "  $SYMLINK_TBI → test_ce.vcf.gz.tbi"

# --- Cleanup ---
rm -f "$TEMP_VCF" "$STRAIN_LIST_FILE" "$HEADER_FILE"

# --- Verify ---
SAMPLE_COUNT=$(bcftools query -l "$OUTPUT_VCF" | wc -l | tr -d ' ')
VARIANT_COUNT=$(bcftools index -n "$OUTPUT_VCF")
FILE_SIZE=$(ls -lh "$OUTPUT_VCF" | awk '{print $5}')
CONTIGS=$(bcftools view -h "$OUTPUT_VCF" | grep '##contig' | sed 's/.*ID=//;s/,.*//')

echo ""
echo "=== C. elegans test VCF generated ==="
echo "  File:       $OUTPUT_VCF"
echo "  Size:       $FILE_SIZE"
echo "  Samples:    $SAMPLE_COUNT"
echo "  Variants:   $VARIANT_COUNT"
echo "  Contigs:    $(echo $CONTIGS | tr '\n' ' ')"
echo ""

# --- Update test_strains_three_species.txt ---
THREE_SPECIES_FILE="${SCRIPT_DIR}/test_strains_three_species.txt"
STRAINS_CSV=$(bcftools query -l "$OUTPUT_VCF" | tr '\n' ',' | sed 's/,$//')
POPA_CSV=$(bcftools query -l "$OUTPUT_VCF" | head -${POP_SPLIT_COUNT} | tr '\n' ',' | sed 's/,$//')
POPB_CSV=$(bcftools query -l "$OUTPUT_VCF" | tail -${POP_SPLIT_COUNT} | tr '\n' ',' | sed 's/,$//')

{
    printf "group\tspecies\tvcf\tms_maf\tms_ld\tstrains\n"
    # Preserve non-ce rows from existing strainfile
    if [[ -f "$THREE_SPECIES_FILE" ]]; then
        awk -F'\t' 'NR>1 && $1 !~ /^ce\./' "$THREE_SPECIES_FILE"
    fi
    printf "ce.test\tc_elegans\tdata/test/test_ce.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$DATE" "$STRAINS_CSV"
    printf "ce.test.popA\tc_elegans\tdata/test/test_ce.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$DATE" "$POPA_CSV"
    printf "ce.test.popB\tc_elegans\tdata/test/test_ce.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$DATE" "$POPB_CSV"
} > "${THREE_SPECIES_FILE}.tmp"
mv "${THREE_SPECIES_FILE}.tmp" "$THREE_SPECIES_FILE"
echo "Updated $THREE_SPECIES_FILE (ce.test, ce.test.popA, ce.test.popB)"

echo "Remember to:"
echo "  1. git add data/test/test_ce.${DATE}.vcf.gz data/test/test_ce.${DATE}.vcf.gz.tbi"
echo "  2. git add data/test/test_strains_three_species.txt"
echo "  3. Update CE_VCF_RELEASE in data/test/release_ids.sh to ${DATE}"
echo ""
echo "Done."
