#!/usr/bin/env bash
# compare_molly.sh — Cross-validate remapping results against the MNEc670k Molly dataset.
#
# This is an optional, standalone validation script specific to the Equine80select array.
# It compares the coordinates and alleles produced by this pipeline against those from
# the MNEc670k.unique_remap.FINAL.csv file (published by McCue et al.).
#
# Usage:
#   bash scripts/compare_molly.sh \
#       -b <remapped_bim_file> \
#       -m <MNEc670k.unique_remap.FINAL.csv> \
#       -o <output_dir>
#
# Inputs:
#   -b  PLINK BIM file from this pipeline (e.g. results/Equine80select_remapped_equCab3.bim)
#   -m  MNEc670k Molly remap CSV (not included in this repo; obtain from original publication)
#   -o  Output directory for comparison results (default: ./molly_comparison)
#
# Outputs (in output_dir/):
#   MNEc670k_remap.tab              Cleaned tab-separated Molly table
#   Equine80select2_molly.bim       EquCab2-source markers remapped via Molly
#   Equine80select3_molly.bim       EquCab3-source markers remapped via Molly
#   Equine80select_molly.bim        Combined Molly BIM
#   molly_remapped_common_snps.txt  SNP IDs in both datasets
#   molly_remapped_matching_snps.txt SNP IDs with identical chr+pos in both
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BIM_FILE=""
MOLLY_FILE=""
OUT_DIR="./molly_comparison"

usage() {
    grep "^#" "$0" | grep -v "^#!" | sed 's/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bim)     BIM_FILE="$2";  shift 2 ;;
        -m|--molly)   MOLLY_FILE="$2"; shift 2 ;;
        -o|--output)  OUT_DIR="$2";   shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$BIM_FILE"   ]] && { echo "ERROR: -b/--bim is required.";   exit 1; }
[[ -z "$MOLLY_FILE" ]] && { echo "ERROR: -m/--molly is required."; exit 1; }
[[ -f "$BIM_FILE"   ]] || { echo "ERROR: BIM file not found: $BIM_FILE";     exit 1; }
[[ -f "$MOLLY_FILE" ]] || { echo "ERROR: Molly file not found: $MOLLY_FILE"; exit 1; }

mkdir -p "$OUT_DIR"

echo "[molly] Cleaning Molly remap table..."
MOLLY_TAB="$OUT_DIR/MNEc670k_remap.tab"
# Keep records with an EquCab3 position; normalise chr names (remove 'chr' prefix, fix unplaced)
awk -F"," '{if($10)print}' "$MOLLY_FILE" \
    | sed 's/chrUn_ref|/Un_/' \
    | sed 's/\.1|,/v1,/' \
    | sed 's/,chr/,/g' \
    | tr ',' '\t' > "$MOLLY_TAB"
echo "[molly] Molly table: $(wc -l < "$MOLLY_TAB") records"

echo "[molly] Generating EquCab2 and EquCab3 source BIM files..."

# For each source map (EquCab2 and EquCab3 markers from the original manifest),
# determine which Molly coordinate space to use (EquCab2 columns 5+7 or EquCab3 columns 6+8)
# by counting how many positions match.

# Build source maps from the BIM file (assume it has chr/snp/0/pos/ref/alt columns)
# Note: This uses the original manifest EquCab2/3 information stored in the BIM SNP IDs;
#       if the BIM file does not carry source-assembly info, provide two separate BIM files.

for f in "$BIM_FILE"; do
    f2="$OUT_DIR/$(basename "${f%.bim}")_molly.bim"
    equ2=$(comm -12 \
        <(awk '{print $1"."$4}' "$f" | sort) \
        <(awk '{print $5"."$7}' "$MOLLY_TAB" | sort) | wc -l)
    equ3=$(comm -12 \
        <(awk '{print $1"."$4}' "$f" | sort) \
        <(awk '{print $6"."$8}' "$MOLLY_TAB" | sort) | wc -l)

    echo "[molly] $(basename "$f"): EquCab2 matches=$equ2, EquCab3 matches=$equ3"

    if [ "$equ2" -gt "$equ3" ]; then
        echo "[molly] → Treating as EquCab2-source; remapping via Molly to EquCab3"
        awk 'BEGIN{FS=OFS="\t"}
            FNR==NR{if($10){a[$5"_"$7]=$6; b[$5"_"$7]=$8; c[$5"_"$7]=$9 FS $10;}next}
            {if(a[$1"_"$4]) print a[$1"_"$4] FS $2 FS $3 FS b[$1"_"$4] FS c[$1"_"$4];}' \
            "$MOLLY_TAB" "$f" | sort -k1,1d -k4,4n > "$f2"
    else
        echo "[molly] → Treating as EquCab3-source; confirming via Molly"
        awk 'BEGIN{FS=OFS="\t"}
            FNR==NR{if($10){a[$6"_"$8]=$6; b[$6"_"$8]=$8; c[$6"_"$8]=$9 FS $10;}next}
            {if(a[$1"_"$4]) print a[$1"_"$4] FS $2 FS $3 FS b[$1"_"$4] FS c[$1"_"$4];}' \
            "$MOLLY_TAB" "$f" | sort -k1,1d -k4,4n > "$f2"
    fi
    echo "[molly] Written: $f2 ($(wc -l < "$f2") records)"
done

echo ""
echo "[molly] Cross-comparison..."
PIPELINE_BIM="$BIM_FILE"
MOLLY_BIM="$OUT_DIR/$(basename "${BIM_FILE%.bim}")_molly.bim"

COMMON="$OUT_DIR/molly_remapped_common_snps.txt"
MATCHING="$OUT_DIR/molly_remapped_matching_snps.txt"

comm -12 \
    <(cut -d$'\t' -f2 "$MOLLY_BIM"   | sort) \
    <(cut -d$'\t' -f2 "$PIPELINE_BIM" | sort) > "$COMMON"

echo "[molly] SNPs in common: $(wc -l < "$COMMON")"

comm -12 \
    <(grep -Fwf "$COMMON" "$MOLLY_BIM"    | sort) \
    <(grep -Fwf "$COMMON" "$PIPELINE_BIM" | sort) | cut -d$'\t' -f2 > "$MATCHING"

echo "[molly] SNPs with identical chr+pos: $(wc -l < "$MATCHING")"

echo ""
echo "[molly] Discrepant markers (first 10):"
paste \
    <(grep -vFwf "$MATCHING" "$COMMON" | grep -Fwf - "$MOLLY_BIM"    | sort -k2,2) \
    <(grep -vFwf "$MATCHING" "$COMMON" | grep -Fwf - "$PIPELINE_BIM" | sort -k2,2) \
    | head -10

echo ""
echo "[molly] Done. Results in: $OUT_DIR"
