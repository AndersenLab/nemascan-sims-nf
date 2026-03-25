#!/usr/bin/env bash
#
# Generate the C. elegans test VCF for nemascan-sims-nf integration testing.
#
# Subsets a CaeNDR C. elegans isotype VCF to the first DEFAULT_STRAIN_COUNT
# strains in its header and chromosomes I, II, V. Removes monomorphic sites,
# strips contig headers for absent chromosomes (III, IV, X), creates a dated
# symlink, and rewrites test_strains.txt so the strain list always matches
# the actual output VCF.
#
# Prerequisites:
#   - bcftools (>= 1.16)
#   - tabix (from htslib)
#   - Source VCF: WI.YYYYMMDD.hard-filter.isotype.vcf.gz
#     Download from CaeNDR: https://elegansvariation.org/data/release/latest
#     See data/test/release_ids.sh for the canonical release date in use.
#
# Usage:
#   ./generate_test_vcf.sh /path/to/WI.YYYYMMDD.hard-filter.isotype.vcf.gz
#
# Output:
#   data/test/test_ce.vcf.gz              (BGZF-compressed VCF)
#   data/test/test_ce.vcf.gz.tbi          (tabix index)
#   data/test/test_ce.YYYYMMDD.vcf.gz     (dated symlink → test_ce.vcf.gz)
#   data/test/test_ce.YYYYMMDD.vcf.gz.tbi (dated symlink → test_ce.vcf.gz.tbi)
#   data/test/test_strains.txt            (rewritten to match output VCF strains)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_VCF="${SCRIPT_DIR}/test_ce.vcf.gz"
STRAIN_FILE="${SCRIPT_DIR}/test_strains.txt"
CHROMOSOMES="I,II,V"
DEFAULT_STRAIN_COUNT=200

# --- Parse arguments ---
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <source_vcf>"
    echo ""
    echo "  source_vcf  Path to WI.YYYYMMDD.hard-filter.isotype.vcf.gz (from CaeNDR)"
    echo "              See data/test/release_ids.sh for the canonical release date."
    exit 1
fi

SOURCE_VCF="$1"

if [[ ! -f "$SOURCE_VCF" ]]; then
    echo "Error: Source VCF not found: $SOURCE_VCF"
    exit 1
fi

# --- Extract release date from source VCF filename ---
VCF_BASENAME="$(basename "$SOURCE_VCF")"
DATE=$(echo "$VCF_BASENAME" | grep -oE '[0-9]{8}' | head -1)
if [[ -z "$DATE" ]]; then
    echo "Error: Cannot extract 8-digit release date from filename: $VCF_BASENAME"
    echo "Rename the file to include a YYYYMMDD date (e.g. WI.20250625.hard-filter.isotype.vcf.gz)"
    exit 1
fi
echo "Detected release date: $DATE"

# --- Check dependencies ---
for cmd in bcftools tabix; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: $cmd not found on PATH"
        exit 1
    fi
done

# --- Build strain list from source VCF header ---
echo "Selecting first ${DEFAULT_STRAIN_COUNT} strains from source VCF header..."
STRAIN_LIST_FILE=$(mktemp)
bcftools query -l "$SOURCE_VCF" | head -${DEFAULT_STRAIN_COUNT} > "$STRAIN_LIST_FILE"
STRAIN_COUNT=$(wc -l < "$STRAIN_LIST_FILE" | tr -d ' ')
echo "Subsetting ${STRAIN_COUNT} strains to chromosomes ${CHROMOSOMES}..."

# --- Subset VCF ---
# Use -S (file) instead of -s (inline) to avoid shell argument length issues
# --min-ac 1 removes monomorphic sites (alleles private to excluded strains)
echo "Running bcftools view (this may take 5-10 minutes for a 7.8 GB source)..."
TEMP_VCF=$(mktemp "${SCRIPT_DIR}/test_XXXXXX")
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

rm -f "$SYMLINK_VCF" "$SYMLINK_TBI"

ln -s test_ce.vcf.gz     "$SYMLINK_VCF"
ln -s test_ce.vcf.gz.tbi "$SYMLINK_TBI"
echo "Created symlinks:"
echo "  $SYMLINK_VCF → test_ce.vcf.gz"
echo "  $SYMLINK_TBI → test_ce.vcf.gz.tbi"

# --- Rewrite test_strains.txt from output VCF ---
# Derive the group name from the existing strainfile if present, else use default.
GROUP="ce.test.200strains"
if [[ -f "$STRAIN_FILE" ]]; then
    EXISTING_GROUP=$(awk -F'\t' 'NR==2 {print $1}' "$STRAIN_FILE")
    [[ -n "$EXISTING_GROUP" ]] && GROUP="$EXISTING_GROUP"
fi
STRAINS_CSV=$(bcftools query -l "$OUTPUT_VCF" | tr '\n' ',' | sed 's/,$//')
{
    printf "group\tspecies\tvcf\tms_maf\tms_ld\tstrains\n"
    printf "%s\tc_elegans\tdata/test/test_ce.%s.vcf.gz\t0.05\t0.8\t%s\n" \
        "$GROUP" "$DATE" "$STRAINS_CSV"
} > "$STRAIN_FILE"
echo "Updated $STRAIN_FILE (group: $GROUP, ${STRAIN_COUNT} strains)"

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

if [[ "$SAMPLE_COUNT" -ne "$STRAIN_COUNT" ]]; then
    echo "WARNING: Expected ${STRAIN_COUNT} samples but got ${SAMPLE_COUNT} in output VCF"
    exit 1
fi

echo "Remember to:"
echo "  1. git add data/test/test_ce.${DATE}.vcf.gz data/test/test_ce.${DATE}.vcf.gz.tbi"
echo "  2. git add data/test/test_strains.txt"
echo "  3. Update CE_VCF_RELEASE in data/test/release_ids.sh to ${DATE}"
echo ""
echo "Done."
