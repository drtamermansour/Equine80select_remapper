#!/usr/bin/env bash
# run_pipeline.sh — Illumina array manifest remapping pipeline.
#
# Remaps markers in an Illumina manifest from their original assembly to a new
# reference genome using a dual-alignment strategy (probe + TopGenomicSeq).
#
# Prerequisites:
#   conda activate remap   (run install.sh first to create this environment)
#
# Usage:
#   bash run_pipeline.sh -i <manifest.csv> -r <reference.fa> -a <assembly> [options]
#
# Required:
#   -i / --manifest      Path to the Illumina manifest CSV
#   -r / --reference     Path to the target reference genome FASTA
#   -a / --assembly      Assembly name for output labels (e.g. equCab3 →
#                        columns like Chr_equCab3); must be explicitly provided
#
# Optional — I/O:
#   -o / --output-dir       Output directory (default: ./output)
#   -t / --threads          Threads for minimap2 (default: 4)
#
# Optional — Filter strictness:
#   --min-anchor             Minimum anchor evidence: dual / topseq (default) / probe
#   --tie-policy             Tie resolution accepted: unique / resolved (default) / avoid_scaffolds
#   --min-refalt-confidence  RefAlt confidence: high / moderate (default) / low
#
# Optional — Thresholds (use 'off' to disable):
#   --min-mapq-topseq N|off  Min MAPQ for TopGenomicSeq alignments (default: 30)
#   --min-mapq-probe  N|off  Min MAPQ for probe alignments (default: off)
#   --max-coord-delta N|off  Remove markers where |probe_coord − CIGAR_coord| > N (default: off)
#
# Optional — Include/exclude (disabled by default):
#   --include-indels          Include indel markers
#   --include-polymorphic     Include markers at polymorphic positions
#   --include-ambiguous-snps  Include ambiguous (A/T, C/G) SNPs
#
# Optional — Operational:
#   --preset      Tune strictness+threshold+include flags together:
#                 strict / default / permissive. Individual flags override.
#   --keep-temp   Keep intermediate FASTA/SAM files
#   --resume      Reuse existing minimap2 SAM files in --temp-dir; skip alignment
#                 only, but still re-run the downstream remap_manifest.py logic
#                 (coordinate resolution, Ref/Alt determination) so the remapped
#                 CSV and remapping_Report.txt reflect the current code. Same
#                 semantics as `remap_manifest.py --resume`.
#   -h / --help   Show this help message
#
# Example:
#   bash run_pipeline.sh \
#       -i manifests/Equine80select_24_20067593_B1.csv \
#       -r equCab3/equCab3_genome.fa \
#       -a equCab3 \
#       -o results/ \
#       -t 8
#
# For HPC/SLURM, see submit_slurm.sh.
#
# Outputs (in output-dir/):
#   remapping/
#     {prefix}_remapped_{assembly}.csv            Full remapped manifest
#     remapping_Report.txt                        Alignment and resolution summary
#     ambiguous_markers.csv                       Markers with ambiguous mapping
#   qc/
#     {prefix}_allele_map_{assembly}.tsv          Allele-map: manifest<->genome crosswalk (main output)
#     {prefix}_remapped_{assembly}.bim            PLINK BIM format
#     {prefix}_remapped_{assembly}.vcf            Final filtered VCF
#     {prefix}_remapped_{assembly}_traced.csv     Full input with per-marker WhyFiltered column
#     QC_Report.txt                               Marker counts at each filter stage
#     diagnostics/                                MAPQ histograms and benchmarks
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ─────────────────────────────────────────────────────────────────
MANIFEST=""
REFERENCE=""
ASSEMBLY=""
OUTPUT_DIR="./output"
THREADS=4
# QC-filter flags are pass-through: empty string = use qc_filter.py's default.
MIN_ANCHOR=""
TIE_POLICY=""
MIN_REFALT_CONFIDENCE=""
MIN_MAPQ_TOPSEQ=""
MIN_MAPQ_PROBE=""
MAX_COORD_DELTA=""
INCLUDE_INDELS=""
INCLUDE_POLYMORPHIC=""
INCLUDE_AMBIGUOUS_SNPS=""
PRESET=""
KEEP_TEMP=""
RESUME=""

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep "^#" "$0" | grep -v "^#!" | sed 's/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--manifest)              MANIFEST="$2";              shift 2 ;;
        -r|--reference)             REFERENCE="$2";             shift 2 ;;
        -a|--assembly)              ASSEMBLY="$2";              shift 2 ;;
        -o|--output-dir)            OUTPUT_DIR="$2";            shift 2 ;;
        -t|--threads)               THREADS="$2";               shift 2 ;;
        --min-anchor)               MIN_ANCHOR="$2";            shift 2 ;;
        --tie-policy)               TIE_POLICY="$2";            shift 2 ;;
        --min-refalt-confidence)    MIN_REFALT_CONFIDENCE="$2"; shift 2 ;;
        --min-mapq-topseq)          MIN_MAPQ_TOPSEQ="$2";       shift 2 ;;
        --min-mapq-probe)           MIN_MAPQ_PROBE="$2";        shift 2 ;;
        --max-coord-delta)          MAX_COORD_DELTA="$2";       shift 2 ;;
        --include-indels)           INCLUDE_INDELS="--include-indels";                 shift ;;
        --include-polymorphic)      INCLUDE_POLYMORPHIC="--include-polymorphic";       shift ;;
        --include-ambiguous-snps)   INCLUDE_AMBIGUOUS_SNPS="--include-ambiguous-snps"; shift ;;
        --preset)                   PRESET="$2";                shift 2 ;;
        --keep-temp)                KEEP_TEMP="--keep-temp";    shift ;;
        --resume)                   RESUME=1;                   shift ;;
        -h|--help)                  usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────────
if [[ -z "$MANIFEST" || -z "$REFERENCE" || -z "$ASSEMBLY" ]]; then
    echo "ERROR: -i/--manifest, -r/--reference, and -a/--assembly are required."
    usage
fi
[[ -f "$MANIFEST"  ]] || { echo "ERROR: Manifest not found: $MANIFEST";   exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: Reference not found: $REFERENCE"; exit 1; }

# Derive prefix from manifest filename
PREFIX="$(basename "$MANIFEST")"
PREFIX="${PREFIX%.csv}"

mkdir -p "$OUTPUT_DIR"
TEMP_DIR="$OUTPUT_DIR/temp"
REMAPPING_DIR="$OUTPUT_DIR/remapping"
QC_DIR="$OUTPUT_DIR/qc"
mkdir -p "$TEMP_DIR" "$REMAPPING_DIR" "$QC_DIR"

REMAPPED_CSV="$REMAPPING_DIR/${PREFIX}_remapped_${ASSEMBLY}.csv"

echo "========================================================"
echo " Manifest Remapping Pipeline"
echo "========================================================"
echo " Manifest:    $MANIFEST"
echo " Reference:   $REFERENCE"
echo " Assembly:    $ASSEMBLY"
echo " Output dir:  $OUTPUT_DIR"
echo " Threads:     $THREADS"
if [[ -n "$PRESET" ]]; then
    echo " Preset:      $PRESET"
fi
echo "========================================================"

# ── Step 1: Index reference if needed ────────────────────────────────────────
echo ""
echo "[pipeline] Step 1: Reference preparation..."
REF_FAI="${REFERENCE}.fai"
if [[ ! -f "$REF_FAI" ]]; then
    echo "[pipeline] Indexing reference with samtools faidx..."
    samtools faidx "$REFERENCE"
fi

# Generate VCF contig definitions
VCF_CONTIGS="$(dirname "$REFERENCE")/vcf_contigs.txt"
if [[ ! -f "$VCF_CONTIGS" ]]; then
    echo "[pipeline] Generating vcf_contigs.txt..."
    awk -F">" 'BEGIN{ref="";reflen=0}
        /^>/{if(ref!="")print "##contig=<ID="ref",length="reflen">";ref=$2;reflen=0;next}
        {reflen+=length($0)}
        END{if(ref!="")print "##contig=<ID="ref",length="reflen">"}' \
        "$REFERENCE" > "$VCF_CONTIGS"
fi

# ── Step 2: Core remapping (Python) ──────────────────────────────────────────
# --resume here mirrors remap_manifest.py's --resume semantics: if SAMs already
# exist in --temp-dir, skip the minimap2 alignment step but always re-run the
# downstream processing (valid-triple filtering, coordinate resolution, Ref/Alt
# determination) so the remapped CSV, remapping_Report.txt, and sidecar CSVs
# reflect the current code. Earlier this branch skipped remap_manifest.py
# entirely when the CSV existed; that made report-text changes invisible after
# re-runs, which was a silent footgun.
echo ""
echo "[pipeline] Step 2: Running remap_manifest.py..."
RESUME_ARG=()
[[ -n "$RESUME" ]] && RESUME_ARG=(--resume)
python "$SCRIPT_DIR/scripts/remap_manifest.py" \
    -i  "$MANIFEST" \
    -r  "$REFERENCE" \
    -o  "$REMAPPED_CSV" \
    -a  "$ASSEMBLY" \
    --threads "$THREADS" \
    --temp-dir "$TEMP_DIR" \
    "${RESUME_ARG[@]}"

# ── Step 3: QC filtering and output generation (Python) ──────────────────────
# Only pass strictness/threshold flags the user actually supplied; otherwise let
# qc_filter.py use its own defaults (or the preset's values).
QC_ARGS=()
[[ -n "$PRESET"                ]] && QC_ARGS+=(--preset "$PRESET")
[[ -n "$MIN_ANCHOR"            ]] && QC_ARGS+=(--min-anchor "$MIN_ANCHOR")
[[ -n "$TIE_POLICY"            ]] && QC_ARGS+=(--tie-policy "$TIE_POLICY")
[[ -n "$MIN_REFALT_CONFIDENCE" ]] && QC_ARGS+=(--min-refalt-confidence "$MIN_REFALT_CONFIDENCE")
[[ -n "$MIN_MAPQ_TOPSEQ"       ]] && QC_ARGS+=(--min-mapq-topseq "$MIN_MAPQ_TOPSEQ")
[[ -n "$MIN_MAPQ_PROBE"        ]] && QC_ARGS+=(--min-mapq-probe  "$MIN_MAPQ_PROBE")
[[ -n "$MAX_COORD_DELTA"       ]] && QC_ARGS+=(--max-coord-delta "$MAX_COORD_DELTA")
[[ -n "$INCLUDE_INDELS"        ]] && QC_ARGS+=("$INCLUDE_INDELS")
[[ -n "$INCLUDE_POLYMORPHIC"   ]] && QC_ARGS+=("$INCLUDE_POLYMORPHIC")
[[ -n "$INCLUDE_AMBIGUOUS_SNPS" ]] && QC_ARGS+=("$INCLUDE_AMBIGUOUS_SNPS")

echo ""
echo "[pipeline] Step 3: Running qc_filter.py..."
python "$SCRIPT_DIR/scripts/qc_filter.py" \
    -i  "$REMAPPED_CSV" \
    -r  "$REFERENCE" \
    -v  "$VCF_CONTIGS" \
    -a  "$ASSEMBLY" \
    -o  "$QC_DIR" \
    --temp-dir "$TEMP_DIR" \
    --prefix   "$PREFIX" \
    "${QC_ARGS[@]}"

# ── Cleanup temp files ────────────────────────────────────────────────────────
if [[ -z "$KEEP_TEMP" ]]; then
    echo "[pipeline] Cleaning up temp files..."
    rm -f "$TEMP_DIR/temp_topseq.fasta" \
          "$TEMP_DIR/temp_probes.fasta" \
          "$TEMP_DIR/temp_topseq.sam" \
          "$TEMP_DIR/temp_probe.sam"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo " Pipeline complete."
echo " Main output: $QC_DIR/${PREFIX}_allele_map_${ASSEMBLY}.tsv"
echo " QC report:   $QC_DIR/QC_Report.txt"
echo " Remap CSV:   $REMAPPED_CSV"
echo "========================================================"
