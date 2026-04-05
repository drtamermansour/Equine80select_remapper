# Scaffold Haplotype Classification Thresholds

Recommended thresholds for classifying unplaced scaffolds (output of
`scripts/scaffold_haplotype_analyzer.py`) as alternative haplotypes of placed chromosomes.

---

## Tier 1 — High confidence alt haplotypes

Safe to treat as chromosomal. All five conditions must be met.

| Metric | Threshold | Rationale |
|---|---|---|
| `identity_pct` | ≥ 99% | Alt haplotypes in a high-quality assembly like EquCab3 diverge <1% from the primary sequence |
| `query_coverage_pct` | ≥ 80% | Most of the scaffold must be accounted for by the alignment |
| `span_to_scaffold_ratio` | ≤ 1.5 | Scaffold maps to a chromosomal window roughly its own size; small SVs push this slightly above 1.0; values >>1 indicate dispersed sequence |
| `max_mapq` | ≥ 40 | Unique placement on one chromosome; MAPQ=0 indicates the scaffold maps equally well to multiple loci (repeat) |
| `n_alignment_blocks` | ≤ 5 | Few blocks = one contiguous region with minor indels; many blocks = fragmented/repetitive alignment |

---

## Tier 2 — Probable alt haplotypes

Worth investigating before use in downstream analysis. All four conditions must be met.

| Metric | Threshold |
|---|---|
| `identity_pct` | ≥ 95% |
| `query_coverage_pct` | ≥ 60% |
| `span_to_scaffold_ratio` | ≤ 3.0 |
| `max_mapq` | ≥ 20 |

---

## Likely repeats / exclude

Scaffolds matching any of these patterns are likely repetitive elements or dispersed duplications,
not alt haplotypes:

- `span_to_scaffold_ratio` > 5 — sequence is dispersed across a large chromosomal window
- `max_mapq` = 0 — equally good alignments exist in multiple chromosomal locations
- `n_alignment_blocks` > 10 — highly fragmented alignment

---

## EquCab3 context

From the EquCab3 run (4668 unplaced scaffolds):

| Category | Count |
|---|---|
| Total scaffolds with ≥1 chromosome hit | 2640 |
| identity ≥ 99% | 722 |
| identity 95–99% | 1497 |
| span_to_scaffold_ratio ≤ 1.5 | 2091 |
| max_mapq = 60 | 1453 |
| Expected Tier 1 | ~500–700 |

---

## Filter script

Use `scripts/filter_scaffold_haplotypes.py` to apply these thresholds programmatically.
Tier 1 defaults are pre-set; all thresholds are adjustable via CLI flags.
