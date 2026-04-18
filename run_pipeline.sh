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
# Optional:
#   -o / --output-dir    Output directory (default: ./output)
#   -t / --threads       Threads for minimap2 (default: 4)
#   --mapq-topseq        Min MAPQ for TopGenomicSeq alignments (default: 30)
#   --mapq-probe         Min MAPQ for probe alignments when >0 (default: 0 = disabled)
#   --coord-delta        Remove markers where |probe_coord − CIGAR_coord| > N (default: -1 = disabled)
#   --coordinate-role    Anchor threshold: High / Moderate (default) / Low
#   --tie-label          Tie threshold: unique / resolved (default) / avoid_scaffolds
#   --refalt-conf        RefAlt confidence: High / Moderate (default) / Low
#   --keep-indels        Include indel markers in outputs (default: excluded)
#   --keep-polymorphic   Keep polymorphic positions (default: excluded)
#   --keep-ambiguous     Keep ambiguous (A/T, C/G) SNPs in outputs (default: excluded)
#   --keep-temp          Keep intermediate FASTA/SAM files
#   --resume             Skip step 2 if remapped CSV and SAM files already exist
#   -h / --help          Show this help message
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
MAPQ_TOPSEQ=30
MAPQ_PROBE=0
COORD_DELTA=-1
COORDINATE_ROLE="Moderate"
TIE_LABEL="resolved"
REFALT_CONF="Moderate"
KEEP_INDELS=""
KEEP_POLYMORPHIC=""
KEEP_AMBIGUOUS=""
KEEP_TEMP=""
RESUME=""

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep "^#" "$0" | grep -v "^#!" | sed 's/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--manifest)      MANIFEST="$2";    shift 2 ;;
        -r|--reference)     REFERENCE="$2";   shift 2 ;;
        -a|--assembly)      ASSEMBLY="$2";    shift 2 ;;
        -o|--output-dir)    OUTPUT_DIR="$2";  shift 2 ;;
        -t|--threads)       THREADS="$2";     shift 2 ;;
        --mapq-topseq)      MAPQ_TOPSEQ="$2"; shift 2 ;;
        --mapq-probe)       MAPQ_PROBE="$2";  shift 2 ;;
        --coord-delta)      COORD_DELTA="$2"; shift 2 ;;
        --coordinate-role)  COORDINATE_ROLE="$2"; shift 2 ;;
        --tie-label)        TIE_LABEL="$2"; shift 2 ;;
        --refalt-conf)      REFALT_CONF="$2"; shift 2 ;;
        --keep-indels)      KEEP_INDELS="--keep-indels"; shift ;;
        --keep-polymorphic) KEEP_POLYMORPHIC="--keep-polymorphic"; shift ;;
        --keep-ambiguous)   KEEP_AMBIGUOUS="--keep-ambiguous"; shift ;;
        --keep-temp)        KEEP_TEMP="--keep-temp"; shift ;;
        --resume)           RESUME=1; shift ;;
        -h|--help)          usage ;;
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
echo " MAPQ TopSeq: $MAPQ_TOPSEQ"
echo " MAPQ Probe:  $MAPQ_PROBE (0 = disabled)"
echo " CoordDelta:  $COORD_DELTA (-1 = disabled)"
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
echo ""
if [[ -n "$RESUME" && -f "$REMAPPED_CSV" ]]; then
    echo "[pipeline] Step 2: Skipping remap_manifest.py (--resume: $REMAPPED_CSV already exists)"
else
    echo "[pipeline] Step 2: Running remap_manifest.py..."
    python "$SCRIPT_DIR/scripts/remap_manifest.py" \
        -i  "$MANIFEST" \
        -r  "$REFERENCE" \
        -o  "$REMAPPED_CSV" \
        -a  "$ASSEMBLY" \
        --threads "$THREADS" \
        --temp-dir "$TEMP_DIR"
fi

# ── Step 3: QC filtering and output generation (Python) ──────────────────────
echo ""
echo "[pipeline] Step 3: Running qc_filter.py..."
python "$SCRIPT_DIR/scripts/qc_filter.py" \
    -i  "$REMAPPED_CSV" \
    -r  "$REFERENCE" \
    -v  "$VCF_CONTIGS" \
    -a  "$ASSEMBLY" \
    -o  "$QC_DIR" \
    --coordinate-role "$COORDINATE_ROLE" \
    --tie-label      "$TIE_LABEL" \
    --refalt-conf    "$REFALT_CONF" \
    --mapq-topseq    "$MAPQ_TOPSEQ" \
    --mapq-probe     "$MAPQ_PROBE" \
    --coord-delta    "$COORD_DELTA" \
    --temp-dir       "$TEMP_DIR" \
    --prefix         "$PREFIX" \
    $KEEP_INDELS \
    $KEEP_POLYMORPHIC \
    $KEEP_AMBIGUOUS

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
