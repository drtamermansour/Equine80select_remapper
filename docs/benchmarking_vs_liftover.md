# Benchmark: our tool vs. liftOver / CrossMap

An apples-to-apples comparison of our probe-sequence-alignment remapper
against the two standard coordinate-lifting tools (UCSC `liftOver` and
`CrossMap`). Built-in ground truth comes from the Equine80select array
itself: v1 was designed on equCab2 (for most markers), v2 re-lists the same
markers at equCab3 coordinates, so v2 provides a trusted equCab3 position for
every v1 marker worth scoring.

## What the experiment does

1. Takes the v1 manifest and filters it to the **GenomeBuild = 2** subset
   — the markers that actually need remapping from equCab2 to equCab3.
   (~75,900 of the 81,974 v1 markers.)
2. Feeds that subset through three independent remappers:
   - **Our tool** — `run_pipeline.sh` (re-aligns probe + TopGenomicSeq to
     equCab3 with minimap2; produces chr, pos, strand, alleles).
   - **liftOver** — UCSC's chain-file coordinate lifter.
   - **CrossMap** — Python/C reimplementation of the same idea.
3. Reads the v2 manifest to build a ground-truth `Name → (chr, pos)` map.
4. For each method × marker, classifies the prediction into one of six verdicts.
5. Writes a per-marker TSV and a summary report.

## Verdict taxonomy

| Verdict | Definition |
|---|---|
| `correct` | predicted chromosome and position exactly match ground truth |
| `wrong_pos_le_10bp` | right chromosome, position off by ≤ 10 bp (gap-placement ambiguity around indels / homopolymers) |
| `wrong_pos_le_1kb`  | right chromosome, off by ≤ 1 kb |
| `wrong_pos_gt_1kb`  | right chromosome, off by > 1 kb (structural disagreement) |
| `wrong_chr` | placed on the wrong chromosome |
| `unmapped` | method could not place the marker |

The three `wrong_pos_*` buckets let you see the **magnitude** of errors, not
just a pass/fail count — useful when judging whether a tool's "wrongs" are
clinically benign or genuinely dangerous.

## How to run

Prerequisites:

- The `remap` conda env (with `liftOver`, `CrossMap`, `minimap2`, `bcftools`, `samtools`).
  `bash install.sh` handles all of this.
- An equCab3 reference FASTA (the pipeline's normal `-r` argument).
- An internet connection on first run, so the orchestrator can auto-download
  the UCSC chain file. It's cached at `genomes/chain/equCab2ToEquCab3.over.chain.gz`.

One command:

```bash
conda activate remap
bash scripts/benchmark/run_benchmark_vs_liftover.sh \
    --v1         manifests/Equine80select_24_20067593_B1.csv \
    --v2         manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --reference  genomes/equCab3/equCab3_genome.fa \
    --output-dir results_liftover_benchmark/ \
    --threads    64
```

Add `--resume` to skip our pipeline's minimap2 step when iterating.

## Output layout

```
results_liftover_benchmark/
├── inputs/
│   ├── v1_equCab2_subset.csv    Illumina manifest, GenomeBuild=2 rows only
│   ├── v1_equCab2.bed           chr-prefixed UCSC BED of equCab2 coordinates
│   └── ground_truth.tsv         Name, chr_equCab3, pos_equCab3 (from v2)
├── ours/
│   ├── remapping/               one minimap2 pass — reused by all three presets
│   ├── qc/                      --preset default (written by run_pipeline.sh)
│   ├── qc_strict/               --preset strict   (re-running qc_filter.py)
│   ├── qc_permissive/           --preset permissive
│   └── temp/
├── liftover/
│   ├── lifted.bed               liftOver's remapped BED
│   └── unmapped.bed             liftOver's rejects
├── crossmap/
│   ├── lifted.bed               CrossMap's remapped BED
│   └── unmapped.bed             CrossMap's rejects
└── report/
    ├── three_way.tsv            per-marker verdicts: 3 presets × ours + liftOver + CrossMap
    └── benchmark_summary.txt    5-row head-to-head + sidebar + disagreements
```

## How the three presets are scored

Our tool is evaluated under all three [presets](cli_reference.md) — `strict`,
`default`, `permissive` — in the same run, so the report shows the full
precision/recall trade-off curve without re-running the slow minimap2 step.

A marker is deemed *placed* by our tool for a given preset iff it appears in
that preset's `allele_map_{assembly}.tsv` (i.e. survived all 11 QC stages under
that preset's settings). Any ground-truth marker absent from the allele map
contributes to that preset's `unmapped` count — this captures both the rare
Stage-1 alignment failures and the QC filter's "deliberate rejections".

The expected behaviour, reading the head-to-head table top-to-bottom:

- **`strict`** should show zero wrong placements but the highest `unmapped`
  count — it keeps only the very safest markers.
- **`default`** should sit between the other two.
- **`permissive`** should dominate `liftOver` / `CrossMap` on every column
  (more `correct`, fewer `wrong_*`, fewer `unmapped`), since it applies only
  the always-on Stage 1/2 filters.

## Reading the report

`benchmark_summary.txt` has three sections:

1. **Head-to-head (coordinates only)** — a table with one row per method and
   columns for each verdict. This is the apples-to-apples comparison.
2. **Sidebar — features only our tool provides** — strand flips, indel
   handling, and allele corrections. liftOver and CrossMap don't claim to do
   any of these, so they're listed here rather than counted against the other
   tools in the head-to-head.
3. **Qualitative disagreement sample** — up to 20 markers where our tool
   "wins" (correct while lift/cross are wrong) and up to 20 where our tool
   "loses", with the raw predictions side-by-side. Useful for spot-checking
   motivating examples.

`three_way.tsv` has one row per ground-truth marker with columns:
`Name, truth_chr, truth_pos, ours_chr, ours_pos, verdict_ours, lift_chr,
lift_pos, verdict_liftover, cross_chr, cross_pos, verdict_crossmap`. Filter it
in any spreadsheet or pandas session to drill into specific failure modes.

## Scope and caveats

- The head-to-head scores **only GenomeBuild=2 markers** (equCab2-designed).
  Our tool is assembly-agnostic (it re-aligns sequences) so it can also
  remap the GenomeBuild=3 subset of v1, but liftOver/CrossMap would need a
  different chain file for that — so they'd be scored unfairly. Excluding
  those markers keeps the comparison clean.
- `correct` is an **exact** position match. If you want to tolerate gap
  placement ambiguity (often ±1 bp around indels), add the
  `wrong_pos_le_10bp` count to the correct count.
- Markers missing from v2 (v1 had them, v2 doesn't) are dropped during input
  preparation — they can't be scored, so they're not penalised against any
  tool. The orchestrator prints how many were dropped.
- The chain file is UCSC's standard `equCab2ToEquCab3.over.chain.gz`. liftOver
  and CrossMap therefore produce identical mapped/unmapped counts; any
  disagreement between them in the per-marker TSV would be a tool bug, not a
  disagreement in the underlying model.
