# Quality-Control Filters

After the alignment step has placed each marker on the new reference,
`qc_filter.py` applies an 11-stage filter cascade. Each stage either passes a
marker through or rejects it; counts at every stage are recorded in
`QC_Report.txt`, and a per-marker `WhyFiltered_{assembly}` column in the
trace CSV records the *first* filter that rejected each marker.

## What each filter does

| # | Stage | What it removes | Default behaviour | Override |
|---|---|---|---|---|
| 1 | Failed markers | Markers that couldn't be placed (`Strand` is `N/A`, meaning unmapped or locus_unresolved) | always on | — |
| 2 | Design conflict | SNPs whose Ref allele doesn't match the reference base, and deletions whose Ref sequence isn't in the genome | always on | — |
| 3 | Coordinate role | Markers placed by less-trusted evidence | filter at `Moderate` (allows `topseq_n_probe` + `topseq_only`) | `--coordinate-role High`/`Low` |
| 4 | Tie label | Markers whose locus required certain tie-break steps (e.g. picking a placed chromosome over a scaffold) | filter at `resolved` (rejects `scaffold_resolved` and `locus_unresolved`) | `--tie-label unique`/`avoid_scaffolds` |
| 5 | Ref/Alt confidence | Markers where the Ref/Alt assignment came from the weaker of the two methods | filter at `Moderate` | `--refalt-conf High`/`Low` |
| 6 | TopGenomicSeq MAPQ | Markers whose context-sequence alignment was low-quality | filter at MAPQ ≥ 30 | `--mapq-topseq N` (0–60) |
| 7 | Probe MAPQ | Same idea, but for the probe alignment | disabled by default | `--mapq-probe N` (0–60) |
| 8 | CoordDelta | Markers where the probe and TopSeq disagree about the exact coordinate | disabled by default | `--coord-delta N` |
| 9 | Indels | Insertion/deletion markers (most genotyping pipelines want SNPs only) | excluded by default | `--keep-indels` |
| 10 | Polymorphic positions | Multiple markers landing at the same position with different Ref/Alt assignments | excluded by default | `--keep-polymorphic` |
| 11 | Ambiguous SNPs | SNPs whose alleles are `{A,T}` or `{C,G}` — strand-ambiguous in downstream tools | excluded by default | `--keep-ambiguous` |

## How to choose stringency

The first three knobs (`--coordinate-role`, `--tie-label`, `--refalt-conf`)
control how confident you want to be in each marker that survives:

- **High / strict** — only the most-trusted markers. Best for precision-sensitive
  applications (e.g. fine-mapping).
- **Moderate / default** — the recommended balance. Drops only what is clearly
  weak.
- **Low / permissive** — keep as many markers as possible. Best when downstream
  software has its own filters and coverage matters more than confidence.

The exempt-when-`NaN` pattern: a marker without a probe alignment carries
`MAPQ_Probe = NaN`, and the `--mapq-probe` filter intentionally lets it
through (NaN means "not measured here", not "MAPQ was zero"). The same applies
to `MAPQ_TopGenomicSeq = NaN` for `probe_only` markers.

## Per-marker traceability

The file `qc/{prefix}_remapped_{assembly}_traced.csv` is the full input
manifest with one extra column, `WhyFiltered_{assembly}`. For markers that
survived all filters this column is empty; otherwise it contains the first
stage label (e.g. `stage_6_mapq_topseq`) that removed the marker.

This lets you:

- Build "rescue" lists of markers a particular filter excluded
- Audit whether the filters are doing what you want
- Compare two pipeline runs at the per-marker level

## Measuring whether a filter is doing its job

When `benchmark_compare.py` is run with `--traced`, the report adds a **QC
Filtration Impact** section that shows, for each stage:

- how many markers it removed,
- what fraction of those were genuinely *non-correct* (per the benchmark),
- how passing-set accuracy improves as each stage is applied.

This is the canonical way to tell whether a filter is rejecting noise or
discarding good data — see [docs/benchmarking.md](benchmarking.md).
