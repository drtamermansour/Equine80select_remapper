# Quality-Control Filters

After the alignment step has placed each marker on the new reference,
`qc_filter.py` applies an 11-stage filter cascade. Each stage either passes a
marker through or rejects it; counts at every stage are recorded in
`QC_Report.txt`, and a per-marker `WhyFiltered_{assembly}` column in the
trace CSV records the *first* filter that rejected each marker.

## What each filter does

| # | Stage | What it removes | Default behaviour | Override |
|---|---|---|---|---|
| 1 | Failed markers | Markers that couldn't be placed (`Strand` is `N/A`: unmapped, `locus_unresolved`, or `refalt_unresolved`) | always on | — |
| 2 | Design conflict | SNPs whose Ref allele doesn't match the reference base, and deletions whose Ref sequence isn't in the genome | always on | — |
| 3 | Min-anchor evidence | Markers placed by less-trusted evidence | filter at `topseq` (allows `topseq_n_probe` + `topseq_only`) | `--min-anchor dual`/`probe` |
| 4 | Tie policy | Markers whose locus required certain tie-break steps (e.g. picking a placed chromosome over a scaffold) | filter at `resolved` (rejects `scaffold_resolved` and `locus_unresolved`) | `--tie-policy unique`/`avoid_scaffolds` |
| 5 | Min-refalt-confidence | Markers where the Ref/Alt assignment came from the weaker of the two methods. Tiers: `high` = `NM_match` + `NM_validated`; `moderate` adds `NM_N/A` + `NM_tied`; `low` adds `NM_only` + `NM_unmatch` + `NM_corrected`. | filter at `moderate` | `--min-refalt-confidence high`/`low` |
| 6 | TopGenomicSeq MAPQ | Markers whose context-sequence alignment was low-quality | filter at MAPQ ≥ 30 | `--min-mapq-topseq N\|off` |
| 7 | Probe MAPQ | Same idea, but for the probe alignment | disabled by default | `--min-mapq-probe N\|off` |
| 8 | Coord-delta | Markers where the probe and TopSeq disagree about the exact coordinate | disabled by default | `--max-coord-delta N\|off` |
| 9 | Indels | Insertion/deletion markers (most genotyping pipelines want SNPs only) | excluded by default | `--include-indels` |
| 10 | Polymorphic positions | Multiple markers landing at the same position with different Ref/Alt assignments | excluded by default | `--include-polymorphic` |
| 11 | Ambiguous SNPs | SNPs whose alleles are `{A,T}` or `{C,G}` — strand-ambiguous in downstream tools | excluded by default | `--include-ambiguous-snps` |

## How to choose stringency

The first three knobs (`--min-anchor`, `--tie-policy`,
`--min-refalt-confidence`) control how confident you want to be in each
marker that survives. For a one-line shortcut, `--preset strict`,
`--preset default`, or `--preset permissive` bundle the strictness,
threshold, and include/exclude flags together; individual flags passed
alongside still override.

- **strict** — only the most-trusted markers. Best for precision-sensitive
  applications (e.g. fine-mapping).
- **default** (recommended) — balanced. Drops only what is clearly weak.
- **permissive** — keep as many markers as possible. Best when downstream
  software has its own filters and coverage matters more than confidence.

The exempt-when-`NaN` pattern: a marker without a probe alignment carries
`MAPQ_Probe = NaN`, and the `--min-mapq-probe` filter intentionally lets
it through (NaN means "not measured here", not "MAPQ was zero"). The same
applies to `MAPQ_TopGenomicSeq = NaN` for `probe_only` markers.

## Per-marker traceability

The file `qc/{prefix}_remapped_{assembly}_traced.csv` is the full input
manifest with one extra column, `WhyFiltered_{assembly}`. For markers that
survived all filters this column is empty; otherwise it contains the first
stage label (e.g. `stage_6_mapq_topseq`) that removed the marker.

The same `stage_N_<slug>` identifiers appear in `QC_Report.txt`'s per-stage
rows, so you can grep the trace CSV for the label you saw in the report.

This lets you:

- Build "rescue" lists of markers a particular filter excluded
- Audit whether the filters are doing what you want
- Compare two pipeline runs at the per-marker level

## "Stage skipped" rows in the report

When a threshold is `off` or an include-flag is set, the corresponding stage
doesn't run but still appears in `QC_Report.txt` — for example
`stage_9_indel_excluded skipped (--include-indels set; would have removed
177)`. The `would have removed` count is the number of markers the stage
*would* have dropped had it run, useful for answering "what would `strict`
do here?" without re-running the pipeline.

## Measuring whether a filter is doing its job

When `benchmark_compare.py` is run with `--traced`, the report adds a **QC
Filtration Impact** section that shows, for each stage:

- how many markers it removed,
- what fraction of those were genuinely *non-correct* (per the benchmark),
- how passing-set accuracy improves as each stage is applied.

This is the canonical way to tell whether a filter is rejecting noise or
discarding good data — see [docs/benchmarking.md](benchmarking.md).
