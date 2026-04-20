#!/usr/bin/env bash
# run_benchmark_vs_liftover.sh — end-to-end orchestrator for the liftOver/CrossMap benchmark.
#
# Runs our pipeline, liftOver, and CrossMap on the GenomeBuild=2 subset of v1,
# then scores all three against v2's equCab3 coordinates.
#
# Usage:
#   bash scripts/benchmark/run_benchmark_vs_liftover.sh \
#       --v1         <v1_manifest.csv>      (mixed-assembly; GenomeBuild=2 rows are the benchmark set)
#       --v2         <v2_manifest.csv>      (all equCab3; ground truth)
#       --reference  <equCab3_genome.fa>
#       --output-dir <results_dir>
#       [--chain     <path>]                (default: genomes/chain/equCab2ToEquCab3.over.chain.gz)
#       [--noAlt-reference <path>]          (FASTA with alt-haplotype scaffolds removed;
#                                            enables the 'permissive_noAlt' arm. If omitted,
#                                            defaults to <REFERENCE_dir>_cleaned/<stem>_no_alt_haplotypes.fa;
#                                            auto-built via the pre-pipeline if that path is missing.)
#       [--threads   N]                     (default: 4)
#       [--resume]                          (skip our pipeline's minimap2 step if already done)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

V1=""
V2=""
REFERENCE=""
OUTPUT_DIR=""
CHAIN=""
THREADS=4
RESUME=""
NOALT_REFERENCE=""

usage() { grep "^#" "$0" | grep -v "^#!" | sed 's/^# \{0,1\}//'; exit 0; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --v1)         V1="$2"; shift 2 ;;
        --v2)         V2="$2"; shift 2 ;;
        --reference)  REFERENCE="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --chain)      CHAIN="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --resume)     RESUME="--resume"; shift ;;
        --noAlt-reference) NOALT_REFERENCE="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$V1" || -z "$V2" || -z "$REFERENCE" || -z "$OUTPUT_DIR" ]] && {
    echo "ERROR: --v1, --v2, --reference, and --output-dir are required."; usage;
}
[[ -f "$V1"        ]] || { echo "ERROR: v1 manifest not found: $V1";    exit 1; }
[[ -f "$V2"        ]] || { echo "ERROR: v2 manifest not found: $V2";    exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: reference not found: $REFERENCE"; exit 1; }

INPUTS_DIR="$OUTPUT_DIR/inputs"
OURS_DIR="$OUTPUT_DIR/ours"
LIFT_DIR="$OUTPUT_DIR/liftover"
CROSS_DIR="$OUTPUT_DIR/crossmap"
REPORT_DIR="$OUTPUT_DIR/report"
mkdir -p "$INPUTS_DIR" "$OURS_DIR" "$LIFT_DIR" "$CROSS_DIR" "$REPORT_DIR"

# ── Chain file ───────────────────────────────────────────────────────────────
if [[ -z "$CHAIN" ]]; then
    CHAIN="$REPO_ROOT/genomes/chain/equCab2ToEquCab3.over.chain.gz"
fi
if [[ ! -f "$CHAIN" ]]; then
    echo "[bench] Chain file not found at $CHAIN — downloading from UCSC..."
    mkdir -p "$(dirname "$CHAIN")"
    wget -q --show-progress -O "$CHAIN" \
        http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/liftOver/equCab2ToEquCab3.over.chain.gz
fi
echo "[bench] Using chain: $CHAIN"

# ── noAlt reference resolution + optional auto-build ─────────────────────────
# The 'permissive_noAlt' arm runs our pipeline against a FASTA with alt-haplotype
# unplaced scaffolds removed (see docs/scaffold_filtering.md). If --noAlt-reference
# is not supplied, default to a sibling '{reference_dir}_cleaned' location whose
# layout matches exclude_alt_haplotypes.py's output. If that path doesn't exist yet,
# run the three-step pre-pipeline to build it from $REFERENCE.
if [[ -z "$NOALT_REFERENCE" ]]; then
    REF_DIR="$(dirname "$REFERENCE")"
    REF_PARENT="$(dirname "$REF_DIR")"
    REF_STEM="$(basename "$REFERENCE")"; REF_STEM="${REF_STEM%.fa}"; REF_STEM="${REF_STEM%.fasta}"
    NOALT_REFERENCE="$REF_PARENT/$(basename "$REF_DIR")_cleaned/${REF_STEM}_no_alt_haplotypes.fa"
fi

if [[ ! -f "$NOALT_REFERENCE" ]]; then
    echo "[bench] noAlt reference not found at $NOALT_REFERENCE — running pre-pipeline to build it."
    PRE_DIR="$(dirname "$NOALT_REFERENCE")/prepipeline"
    mkdir -p "$PRE_DIR" "$(dirname "$NOALT_REFERENCE")"

    python "$REPO_ROOT/scripts/scaffold_haplotype_analyzer.py" \
        -r "$REFERENCE" \
        -o "$PRE_DIR" \
        -t "$THREADS"

    python "$REPO_ROOT/scripts/filter_scaffold_haplotypes.py" \
        -i "$PRE_DIR/scaffold_summary.tsv" \
        -o "$PRE_DIR/alt_haplotype_candidates.tsv"

    python "$REPO_ROOT/scripts/exclude_alt_haplotypes.py" \
        --scaffolds "$PRE_DIR/alt_haplotype_candidates.tsv" \
        --reference "$REFERENCE" \
        --output-dir "$(dirname "$NOALT_REFERENCE")"
fi
[[ -f "$NOALT_REFERENCE" ]] || { echo "ERROR: noAlt reference still missing after auto-build: $NOALT_REFERENCE"; exit 1; }
echo "[bench] Using noAlt reference: $NOALT_REFERENCE"

# ── Step 1: prepare inputs ───────────────────────────────────────────────────
echo ""
echo "[bench] Step 1: prepare benchmark inputs"
python "$SCRIPT_DIR/prepare_benchmark_inputs.py" \
    --v1 "$V1" --v2 "$V2" -o "$INPUTS_DIR"

SUBSET_MANIFEST="$INPUTS_DIR/v1_equCab2_subset.csv"
BED_INPUT="$INPUTS_DIR/v1_equCab2.bed"
GROUND_TRUTH="$INPUTS_DIR/ground_truth.tsv"

# Derive prefix (same rule as run_pipeline.sh).
PREFIX="$(basename "$SUBSET_MANIFEST")"
PREFIX="${PREFIX%.csv}"
OUR_REMAPPED_CSV="$OURS_DIR/remapping/${PREFIX}_remapped_equCab3.csv"
VCF_CONTIGS="$(dirname "$REFERENCE")/vcf_contigs.txt"

# ── Step 2a: our pipeline with --preset default ──────────────────────────────
# This does the reference indexing, vcf_contigs generation, the slow minimap2
# remap, and the QC filter once (with the default preset). The 'default' QC
# output lands in $OURS_DIR/qc/.
echo ""
echo "[bench] Step 2a: run our pipeline (minimap2 + QC, --preset default)"
bash "$REPO_ROOT/run_pipeline.sh" \
    -i "$SUBSET_MANIFEST" \
    -r "$REFERENCE" \
    -a equCab3 \
    -o "$OURS_DIR" \
    -t "$THREADS" \
    --preset default \
    $RESUME

# ── Step 2b/2c: re-run qc_filter.py for the other two presets ────────────────
# We reuse the remapping CSV (the slow step) and just re-apply QC with
# different presets. Each preset's outputs go into their own sibling dir.
for PRESET in strict permissive; do
    TARGET="$OURS_DIR/qc_${PRESET}"
    echo ""
    echo "[bench] Step 2 (preset=$PRESET): re-running qc_filter.py"
    mkdir -p "$TARGET"
    python "$REPO_ROOT/scripts/qc_filter.py" \
        -i  "$OUR_REMAPPED_CSV" \
        -r  "$REFERENCE" \
        -v  "$VCF_CONTIGS" \
        -a  equCab3 \
        -o  "$TARGET" \
        --prefix   "$PREFIX" \
        --temp-dir "$OURS_DIR/temp" \
        --preset   "$PRESET" > "$TARGET/qc_${PRESET}.log" 2>&1
    echo "[bench]   $PRESET  final markers: $(grep '^Final markers' "$TARGET/QC_Report.txt" | awk '{print $NF}')"
done

ALLELE_MAP_DEFAULT="$OURS_DIR/qc/${PREFIX}_allele_map_equCab3.tsv"
ALLELE_MAP_STRICT="$OURS_DIR/qc_strict/${PREFIX}_allele_map_equCab3.tsv"
ALLELE_MAP_PERMISSIVE="$OURS_DIR/qc_permissive/${PREFIX}_allele_map_equCab3.tsv"

# ── Step 2d: our pipeline against the noAlt reference (--preset permissive) ──
# This is a full re-run (new minimap2 alignment) because the reference differs.
# Output lands in a sibling 'ours_noAlt/' so temp/remapping don't collide with
# the raw-genome arms.
OURS_NOALT_DIR="$OUTPUT_DIR/ours_noAlt"
mkdir -p "$OURS_NOALT_DIR"
echo ""
echo "[bench] Step 2d: run our pipeline against noAlt reference (--preset permissive)"
bash "$REPO_ROOT/run_pipeline.sh" \
    -i "$SUBSET_MANIFEST" \
    -r "$NOALT_REFERENCE" \
    -a equCab3 \
    -o "$OURS_NOALT_DIR" \
    -t "$THREADS" \
    --preset permissive \
    $RESUME
ALLELE_MAP_PERMISSIVE_NOALT="$OURS_NOALT_DIR/qc/${PREFIX}_allele_map_equCab3.tsv"

# ── Step 3: liftOver ─────────────────────────────────────────────────────────
echo ""
echo "[bench] Step 3: run liftOver"
liftOver "$BED_INPUT" "$CHAIN" \
    "$LIFT_DIR/lifted.bed" \
    "$LIFT_DIR/unmapped.bed"
echo "[bench] liftOver: $(wc -l < "$LIFT_DIR/lifted.bed") mapped,"\
" $(grep -vc "^#" "$LIFT_DIR/unmapped.bed" || true) unmapped"

# ── Step 4: CrossMap ─────────────────────────────────────────────────────────
echo ""
echo "[bench] Step 4: run CrossMap"
# CrossMap bed writes the mapped BED to the given output path and the unmapped
# BED to <output>.unmap. We normalise the unmapped filename for downstream use.
CrossMap bed "$CHAIN" "$BED_INPUT" "$CROSS_DIR/lifted.bed"
if [[ -f "$CROSS_DIR/lifted.bed.unmap" ]]; then
    mv "$CROSS_DIR/lifted.bed.unmap" "$CROSS_DIR/unmapped.bed"
else
    : > "$CROSS_DIR/unmapped.bed"   # empty file if nothing failed
fi
echo "[bench] CrossMap: $(wc -l < "$CROSS_DIR/lifted.bed") mapped,"\
" $(grep -vc "^#" "$CROSS_DIR/unmapped.bed" || true) unmapped"

# ── Step 5: six-way comparison (4 ours arms + liftOver + CrossMap) ───────────
echo ""
echo "[bench] Step 5: multi-preset comparison"
python "$SCRIPT_DIR/benchmark_vs_liftover.py" \
    --ground-truth      "$GROUND_TRUTH" \
    --ours              "strict:$ALLELE_MAP_STRICT" \
    --ours              "default:$ALLELE_MAP_DEFAULT" \
    --ours              "permissive:$ALLELE_MAP_PERMISSIVE" \
    --ours              "permissive_noAlt:$ALLELE_MAP_PERMISSIVE_NOALT" \
    --liftover-lifted   "$LIFT_DIR/lifted.bed" \
    --liftover-unmapped "$LIFT_DIR/unmapped.bed" \
    --crossmap-lifted   "$CROSS_DIR/lifted.bed" \
    --crossmap-unmapped "$CROSS_DIR/unmapped.bed" \
    -o "$REPORT_DIR"

# ── Step 6: explain every permissive-loss marker against the documented categories
echo ""
echo "[bench] Step 6: explain permissive losses (cross-reference vs docs/why_we_right.md + docs/cant_remap.md)"
python "$SCRIPT_DIR/explain_permissive_losses.py" \
    --three-way   "$REPORT_DIR/three_way.tsv" \
    --v1-manifest "$V1" \
    --remapped    "$OUR_REMAPPED_CSV" \
    --reference   "$REFERENCE" \
    -o            "$REPORT_DIR"

echo ""
echo "[bench] =========================================="
echo "[bench]  Benchmark complete."
echo "[bench]  Summary:    $REPORT_DIR/benchmark_summary.txt"
echo "[bench]  Per-marker: $REPORT_DIR/three_way.tsv"
echo "[bench]  Losses:     $REPORT_DIR/permissive_losses_explained.md"
echo "[bench] =========================================="
