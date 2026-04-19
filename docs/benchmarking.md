# Benchmarking Remapping Accuracy

If you're moving a manifest from assembly A → assembly B and there happens
to exist a separate, independently-curated B-native version of the same
array, you can use the second manifest as ground truth and measure how
accurate the remapping was.

This is exactly what the original Equine80select v1 (EquCab2) →
Equine80select v2 (EquCab3) project did.

## The benchmark scripts

| Script | What it does |
|---|---|
| [`scripts/benchmark_compare.py`](../scripts/benchmark_compare.py) | The main benchmark. Compares the remapped CSV against a ground-truth manifest position-by-position. With `--reference`, also runs an explanatory layer that classifies each non-correct marker. With `--traced`, adds a QC-filtration impact section. |
| [`scripts/benchmark_cigar_vs_probe.py`](../scripts/benchmark_cigar_vs_probe.py) | A diagnostic three-way comparison: probe-CIGAR vs TopSeq-CIGAR vs the final chosen coordinate. Useful for tuning the CoordSource decision rule. |
| [`scripts/benchmark/run_benchmark_vs_liftover.sh`](../scripts/benchmark/run_benchmark_vs_liftover.sh) | **Head-to-head vs. UCSC liftOver and CrossMap** on the equCab2→equCab3 example. See [docs/benchmarking_vs_liftover.md](benchmarking_vs_liftover.md). |

## Quick run

```bash
python scripts/benchmark_compare.py \
    --manifest  ground_truth.csv \
    --remapped  results/remapping/{prefix}_remapped_equCab3.csv \
    --assembly  equCab3 \
    --reference equCab3/equCab3_genome.fa \
    --traced    results/qc/{prefix}_remapped_equCab3_traced.csv \
    --output-dir results/benchmark/
```

The four most important arguments:

- `--manifest` — the ground-truth manifest (must contain `Chr`, `MapInfo`,
  `RefStrand`, `SNP`, `TopGenomicSeq`).
- `--remapped` — the pipeline's output CSV.
- `--reference` (optional) — enables the **explanatory verdict layer** that
  uses sequence context to label each non-correct marker.
- `--traced` (optional) — enables the **QC filtration impact** section that
  shows how each filter stage affected the benchmark outcome distribution.

## What the report contains

Three layers of detail, each opt-in via the corresponding flag:

### 1. Headline (always)

```
correct                 81,943 ( 99.7%)
coord_off                  114 (  0.1%)
wrong_chr                   31 (  0.0%)
unmapped                   134 (  0.2%)
```

Plus a 3-dimension breakdown by `anchor × tie × benchmark result`.

### 2. Explanatory verdict (with `--reference`)

For each non-correct marker, the script extracts ~20 bp flanking context from
the reference and tries to match it against the manifest's `TopGenomicSeq`.
This produces a verdict explaining *who* is wrong:

| Verdict | Meaning |
|---|---|
| `manifest_strand_wrong` | Pipeline placed the marker correctly; manifest's `RefStrand` disagrees with the actual reference orientation. |
| `manifest_coord_wrong` | Pipeline coord differs from manifest's, but context matches at *our* coord — manifest is stale. |
| `pipeline_wrong_locus` | Pipeline placed the marker on the **right chromosome** but context fails at our coord — we picked the wrong position within that chromosome. |
| `pipeline_wrong_chr` | Pipeline placed the marker on the **wrong chromosome**; context can't match by construction. |
| `pipeline_unmapped` | Pipeline assigned no position (`Chr=0` / `Strand=N/A`). Not a "wrong locus" — no locus at all. |
| `pipeline_wrong_strand` | Context matches in *reverse* orientation; we placed the right locus but flipped the strand. |
| `ambiguous_snp` | A/T or C/G SNP — strand is undetermined from alleles alone. |
| `probe_only_inconclusive` | Known-weak class; segregated to keep noise out of the headline. |
| `unresolved` | Manual investigation required. |

Worked examples in [why_we_right.md](why_we_right.md) (manifest wrong) and
[cant_remap.md](cant_remap.md) (genuinely un-remappable markers).

### 3. QC filtration impact (with `--traced`)

Crosses the benchmark verdict with the `WhyFiltered_{assembly}` column from
`qc_filter.py`'s trace CSV, producing:

- **Confusion matrix** — for each filter stage, how many of its removals
  were correct, coord_off, wrong_chr, etc.
- **Passing-set accuracy** — what fraction of markers that survived QC are
  truly correct.
- **Per-stage precision** — `% of removed = non-correct` per stage. A stage
  with 100% precision removes only genuine errors; a stage with low precision
  removes correct markers as well (often a quality-vs-coverage trade-off,
  e.g. ambiguous-SNP exclusion).
- **Cumulative accuracy gain** — passing-set accuracy after each stage is
  applied in order. Tells you which stages have the biggest payoff.
- **False positives** — correct markers QC removed (with the first 10 listed
  for inspection).
- **False negatives** — non-correct markers QC let through (with verdict
  drill-down).

This last layer is the key tool for tuning filter thresholds.

## Three-way coordinate comparison

```bash
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  ground_truth.csv \
    --remapped  results/remapping/{prefix}_remapped_equCab3.csv \
    --assembly  equCab3
```

Reports per-marker accuracy of the probe-CIGAR coordinate, the TopSeq-CIGAR
coordinate, and the final chosen coordinate, stratified by `CoordDelta`
bucket. Used to validate that the `CoordDelta ≥ 2 → use TopSeq CIGAR` rule
in `remap_manifest.py` is empirically correct.

## Output files

All outputs are timestamped:

| File | What's in it |
|---|---|
| `benchmark_{ts}.tsv` | One row per benchmarked marker with every diagnostic signal |
| `benchmark_{ts}_chrY.tsv` | Chr=Y markers (excluded from headline metrics) |
| `benchmark_{ts}_chr0.tsv` | Chr=0 / unplaced markers |
| `benchmark_{ts}_report.txt` | Human-readable report (headline + explanatory + QC impact) |
| `benchmark_{ts}_diff.txt` | Category transitions vs `--baseline` (only when supplied) |
