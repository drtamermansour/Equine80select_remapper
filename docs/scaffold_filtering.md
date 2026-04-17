# Pre-Pipeline: Scaffold Haplotype Filtering

Modern reference assemblies often include hundreds of unplaced scaffolds.
Many of these are *alternative haplotypes* of regions that already appear on
the placed chromosomes — for example, a heterozygous Major Histocompatibility
Complex (MHC) allele assembled separately from the chromosome-resident copy.

When both copies are present in the reference, every probe matching that
region will multi-map between the placed chromosome and its alt-haplotype
scaffold. minimap2 then reports a low MAPQ, and the marker either gets
removed by the MAPQ filter or placed semi-randomly on whichever copy scored
higher.

The fix is to identify and remove these alt-haplotype scaffolds **before**
running the main pipeline.

---

## The three-step workflow

### Step 1 — Characterise unplaced scaffolds

```bash
python scripts/scaffold_haplotype_analyzer.py \
    -r equCab3/equCab3_genome.fa \
    -o scaffold_haplotype_analysis/ \
    --threads 8
```

Aligns every unplaced scaffold to the placed chromosomes with
`minimap2 -x asm5` and writes `scaffold_summary.tsv` with per-scaffold
statistics:

| Column | Meaning |
|---|---|
| `identity_pct` | % identity of the best alignment |
| `query_coverage_pct` | Fraction of the scaffold that aligned |
| `span_to_scaffold_ratio` | Aligned span vs scaffold length (close to 1 = clean alignment, > 5 = chimeric / multi-mapping) |
| `max_mapq` | Best MAPQ across all alignments |
| `n_alignment_blocks` | Number of distinct alignment blocks |

### Step 2 — Pick the alt-haplotype candidates

```bash
python scripts/filter_scaffold_haplotypes.py \
    -i scaffold_haplotype_analysis/scaffold_summary.tsv \
    -o scaffold_haplotype_analysis/alt_haplotype_candidates.tsv
```

Applies a threshold on each statistic. Defaults (Tier 1, high confidence) are:

| Flag | Default | Meaning |
|---|---|---|
| `--min-identity` | 99.0 | Minimum `identity_pct` |
| `--min-query-cov` | 80.0 | Minimum `query_coverage_pct` |
| `--max-span-ratio` | 1.5 | Maximum `span_to_scaffold_ratio` |
| `--min-mapq` | 40 | Minimum `max_mapq` |
| `--max-blocks` | 5 | Maximum `n_alignment_blocks` |

Stricter or looser tiers are documented in
[scaffold_haplotype_thresholds.md](scaffold_haplotype_thresholds.md).

### Step 3 — Build a cleaned reference

```bash
python scripts/exclude_alt_haplotypes.py \
    --scaffolds scaffold_haplotype_analysis/alt_haplotype_candidates.tsv \
    --reference equCab3/equCab3_genome.fa \
    --output-dir equCab3_cleaned/
```

Removes the identified scaffolds from the FASTA and writes:

| File | Description |
|---|---|
| `{stem}_no_alt_haplotypes.fa` | Cleaned FASTA |
| `{stem}_no_alt_haplotypes.fa.fai` | samtools index |
| `exclusion_report.txt` | Counts of excluded vs retained sequences |

Use the cleaned FASTA as the `-r` input to `run_pipeline.sh`.

---

## When this matters

The scaffold filter is **optional**. For most genotyping arrays you'll see
the impact only on a small fraction of markers (typically <0.5%) — those are
the ones documented in [cant_remap.md](cant_remap.md) under the
`both_match_duplicate_locus` category. If you don't run the pre-filter, those
markers will land on whichever copy minimap2 scored highest, which is
arbitrary.

For arrays where every marker counts (fine-mapping panels, low-density
QTL panels), the scaffold filter is worth running.
