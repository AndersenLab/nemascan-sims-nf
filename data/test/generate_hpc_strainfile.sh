#!/usr/bin/env bash
#
# Generate HPC strainfiles for nemascan-sims-nf integration testing on Rockfish.
#
# Reads strain names from full CaeNDR VCFs using bcftools query -l and writes:
#   data/test/hpc_strains.txt                — single CE group (ce.hpc)
#   data/test/hpc_strains_three_species.txt  — 9 groups: ce.hpc + ce.hpc.popA/B
#                                               cb.hpc + cb.hpc.popA/B
#                                               ct.hpc + ct.hpc.popA/B
#
# Prerequisites:
#   - bcftools (>= 1.16) on PATH
#     On Rockfish: module load bcftools && export PATH
#   - Tabix-indexed VCF files (.tbi index must exist alongside .vcf.gz)
#   - Access to the Rockfish /vast/eande106 filesystem
#
# Usage:
#   bash data/test/generate_hpc_strainfile.sh <ce_vcf> [cb_vcf ct_vcf]
#
#   ce_vcf   Required. Absolute path to C. elegans full VCF.
#   cb_vcf   Optional (must pair with ct_vcf). C. briggsae full VCF.
#   ct_vcf   Optional (must pair with cb_vcf). C. tropicalis full VCF.
#
# If only ce_vcf is provided, only hpc_strains.txt is updated.
# If all three are provided, both output files are updated.
# CB and CT must be provided together — a partial three-species file is invalid.
#
# Canonical Rockfish VCF paths (dates match data/test/release_ids.sh):
#   /vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.isotype.vcf.gz
#   /vast/eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz
#   /vast/eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz
#
# After running, commit the updated strainfiles:
#   git add data/test/hpc_strains.txt data/test/hpc_strains_three_species.txt
#   git commit -m "Update HPC strainfiles from Rockfish VCFs"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HPC_STRAINS_FILE="${SCRIPT_DIR}/hpc_strains.txt"
HPC_THREE_SPECIES_FILE="${SCRIPT_DIR}/hpc_strains_three_species.txt"
POP_SPLIT_COUNT=100   # strains per popA / popB (must be <= total strains in each VCF / 2)

# --- Parse arguments ---
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <ce_vcf> [cb_vcf ct_vcf]"
    echo ""
    echo "  ce_vcf   Required. Path to C. elegans WI.YYYYMMDD.hard-filter.isotype.vcf.gz"
    echo "  cb_vcf   Optional (pair with ct_vcf). C. briggsae VCF path."
    echo "  ct_vcf   Optional (pair with cb_vcf). C. tropicalis VCF path."
    echo ""
    echo "Canonical Rockfish paths:"
    echo "  /vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.isotype.vcf.gz"
    echo "  /vast/eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz"
    echo "  /vast/eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz"
    exit 1
fi

CE_VCF="$1"
CB_VCF="${2:-}"
CT_VCF="${3:-}"

# CB and CT must be provided together
if [[ -n "$CB_VCF" ]] && [[ -z "$CT_VCF" ]]; then
    echo "Error: cb_vcf provided without ct_vcf. Supply both or neither." >&2
    exit 1
fi
if [[ -z "$CB_VCF" ]] && [[ -n "$CT_VCF" ]]; then
    echo "Error: ct_vcf provided without cb_vcf. Supply both or neither." >&2
    exit 1
fi

# --- Check dependencies ---
if ! command -v bcftools &>/dev/null; then
    echo "Error: bcftools not found on PATH" >&2
    echo "       On Rockfish: module load bcftools && export PATH" >&2
    exit 1
fi

# --- Helper: validate VCF exists and is tabix-indexed ---
check_vcf() {
    local vcf_path="$1"
    local label="$2"
    if [[ ! -f "$vcf_path" ]]; then
        echo "Error: $label VCF not found: $vcf_path" >&2
        exit 1
    fi
    if [[ ! -f "${vcf_path}.tbi" ]]; then
        echo "Error: $label VCF index not found: ${vcf_path}.tbi" >&2
        echo "       Run: tabix -p vcf \"$vcf_path\"" >&2
        exit 1
    fi
}

# --- Helper: extract 8-digit date from VCF filename ---
extract_date() {
    local vcf_path="$1"
    local basename
    basename="$(basename "$vcf_path")"
    local date
    date=$(echo "$basename" | grep -oE '[0-9]{8}' | head -1)
    if [[ -z "$date" ]]; then
        echo "Error: Cannot extract 8-digit release date from: $basename" >&2
        echo "       Expected filename pattern: WI.YYYYMMDD.hard-filter.isotype.vcf.gz" >&2
        exit 1
    fi
    echo "$date"
}

# --- Process C. elegans (required) ---
check_vcf "$CE_VCF" "C. elegans"
CE_DATE=$(extract_date "$CE_VCF")
echo "C. elegans: $CE_VCF (release: $CE_DATE)"

CE_STRAINS_CSV=$(bcftools query -l "$CE_VCF" | tr '\n' ',' | sed 's/,$//')
CE_STRAIN_COUNT=$(bcftools query -l "$CE_VCF" | wc -l | tr -d ' ')
CE_POPA_CSV=$(bcftools query -l "$CE_VCF" | head -"${POP_SPLIT_COUNT}" | tr '\n' ',' | sed 's/,$//')
CE_POPB_CSV=$(bcftools query -l "$CE_VCF" | tail -"${POP_SPLIT_COUNT}" | tr '\n' ',' | sed 's/,$//')
echo "  Strains: ${CE_STRAIN_COUNT}"

# Write hpc_strains.txt (single CE group)
{
    printf "group\tspecies\tvcf\tms_maf\tms_ld\tstrains\n"
    printf "ce.hpc\tc_elegans\t%s\t0.05\t0.8\t%s\n" "$CE_VCF" "$CE_STRAINS_CSV"
} > "${HPC_STRAINS_FILE}.tmp"
mv "${HPC_STRAINS_FILE}.tmp" "$HPC_STRAINS_FILE"
echo "Written: $HPC_STRAINS_FILE"

# --- If no CB/CT provided, stop here ---
if [[ -z "$CB_VCF" ]]; then
    echo ""
    echo "Only CE VCF provided — hpc_strains_three_species.txt not updated."
    echo "Supply all three VCF paths to also update the three-species strainfile."
    echo ""
    echo "Remember to:"
    echo "  git add data/test/hpc_strains.txt"
    echo "  git commit -m 'Update HPC strainfiles from Rockfish VCFs'"
    echo ""
    echo "Done."
    exit 0
fi

# --- Process C. briggsae ---
check_vcf "$CB_VCF" "C. briggsae"
CB_DATE=$(extract_date "$CB_VCF")
echo "C. briggsae: $CB_VCF (release: $CB_DATE)"

CB_STRAINS_CSV=$(bcftools query -l "$CB_VCF" | tr '\n' ',' | sed 's/,$//')
CB_STRAIN_COUNT=$(bcftools query -l "$CB_VCF" | wc -l | tr -d ' ')
CB_POPA_CSV=$(bcftools query -l "$CB_VCF" | head -"${POP_SPLIT_COUNT}" | tr '\n' ',' | sed 's/,$//')
CB_POPB_CSV=$(bcftools query -l "$CB_VCF" | tail -"${POP_SPLIT_COUNT}" | tr '\n' ',' | sed 's/,$//')
echo "  Strains: ${CB_STRAIN_COUNT}"

# --- Process C. tropicalis ---
check_vcf "$CT_VCF" "C. tropicalis"
CT_DATE=$(extract_date "$CT_VCF")
echo "C. tropicalis: $CT_VCF (release: $CT_DATE)"

CT_STRAINS_CSV=$(bcftools query -l "$CT_VCF" | tr '\n' ',' | sed 's/,$//')
CT_STRAIN_COUNT=$(bcftools query -l "$CT_VCF" | wc -l | tr -d ' ')
CT_POPA_CSV=$(bcftools query -l "$CT_VCF" | head -"${POP_SPLIT_COUNT}" | tr '\n' ',' | sed 's/,$//')
CT_POPB_CSV=$(bcftools query -l "$CT_VCF" | tail -"${POP_SPLIT_COUNT}" | tr '\n' ',' | sed 's/,$//')
echo "  Strains: ${CT_STRAIN_COUNT}"

# Write hpc_strains_three_species.txt (9 groups: full + popA + popB per species)
{
    printf "group\tspecies\tvcf\tms_maf\tms_ld\tstrains\n"
    printf "ce.hpc\tc_elegans\t%s\t0.05\t0.8\t%s\n"       "$CE_VCF" "$CE_STRAINS_CSV"
    printf "ce.hpc.popA\tc_elegans\t%s\t0.05\t0.8\t%s\n"   "$CE_VCF" "$CE_POPA_CSV"
    printf "ce.hpc.popB\tc_elegans\t%s\t0.05\t0.8\t%s\n"   "$CE_VCF" "$CE_POPB_CSV"
    printf "cb.hpc\tc_briggsae\t%s\t0.05\t0.8\t%s\n"       "$CB_VCF" "$CB_STRAINS_CSV"
    printf "cb.hpc.popA\tc_briggsae\t%s\t0.05\t0.8\t%s\n"  "$CB_VCF" "$CB_POPA_CSV"
    printf "cb.hpc.popB\tc_briggsae\t%s\t0.05\t0.8\t%s\n"  "$CB_VCF" "$CB_POPB_CSV"
    printf "ct.hpc\tc_tropicalis\t%s\t0.05\t0.8\t%s\n"     "$CT_VCF" "$CT_STRAINS_CSV"
    printf "ct.hpc.popA\tc_tropicalis\t%s\t0.05\t0.8\t%s\n" "$CT_VCF" "$CT_POPA_CSV"
    printf "ct.hpc.popB\tc_tropicalis\t%s\t0.05\t0.8\t%s\n" "$CT_VCF" "$CT_POPB_CSV"
} > "${HPC_THREE_SPECIES_FILE}.tmp"
mv "${HPC_THREE_SPECIES_FILE}.tmp" "$HPC_THREE_SPECIES_FILE"
echo "Written: $HPC_THREE_SPECIES_FILE"

echo ""
echo "=== Summary ==="
printf "  CE strains: %s (release %s)\n" "$CE_STRAIN_COUNT" "$CE_DATE"
printf "  CB strains: %s (release %s)\n" "$CB_STRAIN_COUNT" "$CB_DATE"
printf "  CT strains: %s (release %s)\n" "$CT_STRAIN_COUNT" "$CT_DATE"
echo ""
echo "Remember to:"
echo "  git add data/test/hpc_strains.txt data/test/hpc_strains_three_species.txt"
echo "  git commit -m 'Update HPC strainfiles from Rockfish VCFs'"
echo ""
echo "Done."
