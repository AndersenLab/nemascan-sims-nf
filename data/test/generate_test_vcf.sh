#!/usr/bin/env bash
#
# Generate the test VCF for nemascan-sims-nf integration testing.
#
# Subsets the CaeNDR isotype reference VCF to the strains listed in
# test_strains.txt and chromosomes I, II, V. Removes monomorphic sites
# and strips contig headers for absent chromosomes (III, IV, X).
#
# Prerequisites:
#   - bcftools (>= 1.16)
#   - tabix (from htslib)
#   - Source VCF: WI.20220216.hard-filter.isotype.vcf.gz
#     Download from CaeNDR: https://elegansvariation.org/data/release/latest
#
# Usage:
#   ./generate_test_vcf.sh /path/to/WI.20220216.hard-filter.isotype.vcf.gz
#
# Output:
#   data/test/test.vcf.gz      (BGZF-compressed VCF)
#   data/test/test.vcf.gz.tbi  (tabix index)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_VCF="${SCRIPT_DIR}/test.vcf.gz"
STRAIN_FILE="${SCRIPT_DIR}/test_strains.txt"
CHROMOSOMES="I,II,V"

# --- Parse arguments ---
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <source_vcf>"
    echo ""
    echo "  source_vcf  Path to WI.20220216.hard-filter.isotype.vcf.gz"
    echo "              (550 strains, 7.8 GB, from CaeNDR)"
    exit 1
fi

SOURCE_VCF="$1"

if [[ ! -f "$SOURCE_VCF" ]]; then
    echo "Error: Source VCF not found: $SOURCE_VCF"
    exit 1
fi

if [[ ! -f "$STRAIN_FILE" ]]; then
    echo "Error: Strain file not found: $STRAIN_FILE"
    exit 1
fi

# --- Check dependencies ---
for cmd in bcftools tabix; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: $cmd not found on PATH"
        exit 1
    fi
done

# --- Extract strain list from test_strains.txt ---
# Format: "{pop_name} {strain1,strain2,...}" — extract the comma-separated part
STRAIN_LIST_FILE=$(mktemp)
awk '{print $2}' "$STRAIN_FILE" | tr ',' '\n' > "$STRAIN_LIST_FILE"
STRAIN_COUNT=$(wc -l < "$STRAIN_LIST_FILE" | tr -d ' ')
echo "Subsetting ${STRAIN_COUNT} strains to chromosomes ${CHROMOSOMES}..."

# --- Subset VCF ---
# Use -S (file) instead of -s (inline) to avoid shell argument length issues
# --min-ac 1 removes monomorphic sites (alleles private to excluded strains)
echo "Running bcftools view (this may take 5-10 minutes for a 7.8 GB source)..."
TEMP_VCF=$(mktemp "${SCRIPT_DIR}/test_XXXXXX.vcf.gz")
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

# --- Cleanup ---
rm -f "$TEMP_VCF" "$STRAIN_LIST_FILE" "$HEADER_FILE"

# --- Verify ---
SAMPLE_COUNT=$(bcftools query -l "$OUTPUT_VCF" | wc -l | tr -d ' ')
VARIANT_COUNT=$(bcftools index -n "$OUTPUT_VCF")
FILE_SIZE=$(ls -lh "$OUTPUT_VCF" | awk '{print $5}')
CONTIGS=$(bcftools view -h "$OUTPUT_VCF" | grep '##contig' | sed 's/.*ID=//;s/,.*//')

echo ""
echo "=== Test VCF generated ==="
echo "  File:       $OUTPUT_VCF"
echo "  Size:       $FILE_SIZE"
echo "  Samples:    $SAMPLE_COUNT"
echo "  Variants:   $VARIANT_COUNT"
echo "  Contigs:    $(echo $CONTIGS | tr '\n' ' ')"
echo ""

if [[ "$SAMPLE_COUNT" -ne "$STRAIN_COUNT" ]]; then
    echo "WARNING: Expected $STRAIN_COUNT samples but got $SAMPLE_COUNT"
    echo "Some strains in $STRAIN_FILE may not be in the source VCF."
    exit 1
fi

echo "Done."
