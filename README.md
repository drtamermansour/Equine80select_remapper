# Array Manifest Remapper

A computational pipeline for remapping Illumina genotyping array manifests between reference genome assemblies. Uses a **dual-alignment strategy** — aligning both the short physical probe (50 bp) and the longer `TopGenomicSeq` context sequence — to ensure high-fidelity coordinate conversion, correct strand assignment, and precise Ref/Alt allele determination consistent with VCF standards.

Originally developed to remap the **Equine80select** array from EquCab2 to EquCab3.

---

## Prerequisites

- **Conda** or **Mamba** package manager
- **Reference genome FASTA** indexed with `samtools faidx`

## Setup

```bash
git clone git@github.com:drtamermansour/Equine80select_remapper.git
cd Equine80select_remapper
bash install.sh
conda activate remap
```

---

## Pre-pipeline: Reference Preparation (Optional)

Modern genome assemblies include unplaced scaffolds that are often alternative haplotypes of placed chromosomes. Including them in the reference causes ambiguous multi-mapping, reducing remapping confidence. Three scripts handle scaffold characterisation and exclusion **before** running the main pipeline.

### Step 1 — Characterise scaffolds

```bash
python scripts/scaffold_haplotype_analyzer.py \
    -r equCab3/equCab3_genome.fa \
    -o scaffold_haplotype_analysis/ \
    --threads 8
```

Aligns all unplaced scaffolds to the placed chromosomes using `minimap2 -x asm5` and produces `scaffold_summary.tsv` with per-scaffold alignment statistics (`identity_pct`, `query_coverage_pct`, `span_to_scaffold_ratio`, `max_mapq`, `n_alignment_blocks`).

### Step 2 — Filter to alt-haplotype candidates

```bash
python scripts/filter_scaffold_haplotypes.py \
    -i scaffold_haplotype_analysis/scaffold_summary.tsv \
    -o scaffold_haplotype_analysis/alt_haplotype_candidates.tsv
```

Applies threshold filters to select scaffolds likely to be alternative haplotypes. Default thresholds are Tier 1 (high confidence). See **[docs/scaffold_haplotype_thresholds.md](docs/scaffold_haplotype_thresholds.md)** for the full threshold rationale and multiple strictness tiers.

Key CLI flags:

| Flag | Default | Meaning |
|---|---|---|
| `--min-identity` | 99.0 | Minimum `identity_pct` |
| `--min-query-cov` | 80.0 | Minimum `query_coverage_pct` |
| `--max-span-ratio` | 1.5 | Maximum `span_to_scaffold_ratio` |
| `--min-mapq` | 40 | Minimum `max_mapq` |
| `--max-blocks` | 5 | Maximum `n_alignment_blocks` |

### Step 3 — Build cleaned reference

```bash
python scripts/exclude_alt_haplotypes.py \
    --scaffolds scaffold_haplotype_analysis/alt_haplotype_candidates.tsv \
    --reference equCab3/equCab3_genome.fa \
    --output-dir equCab3_cleaned/
```

Removes the identified scaffolds from the FASTA and writes an indexed cleaned reference. Outputs:

| File | Description |
|---|---|
| `{stem}_no_alt_haplotypes.fa` | Cleaned FASTA with alt-haplotype scaffolds removed |
| `{stem}_no_alt_haplotypes.fa.fai` | samtools index |
| `exclusion_report.txt` | Count of excluded / retained sequences; lists any scaffold IDs not found in reference |

Use the cleaned FASTA as the `-r` input to the main pipeline.

---

## Running the Pipeline

```bash
bash run_pipeline.sh \
    -i backup_original/Equine80select_24_20067593_B1.csv \
    -r equCab3/equCab3_genome.fa \
    -a equCab3 \
    -o results/
```

### All Options

| Flag | Default | Description |
|---|---|---|
| `-i / --manifest` | *(required)* | Path to the Illumina manifest CSV |
| `-r / --reference` | *(required)* | Path to the target reference genome FASTA |
| `-a / --assembly` | derived from FASTA filename | Assembly name used to label outputs |
| `-o / --output-dir` | `./output` | Output directory |
| `-t / --threads` | `4` | Threads for minimap2 |
| `--mapq-topseq` | `30` | Minimum MAPQ for TopGenomicSeq alignments |
| `--mapq-probe` | `0` (disabled) | Minimum MAPQ for probe alignments |
| `--coord-delta` | `-1` (disabled) | Remove markers where `\|probe_coord − CIGAR_coord\| > N` and all markers where `anchor_{assembly} == "topseq_only"` |
| `--exclude-indels` | off | Remove indel markers from all outputs (VCF, BIM, map file) |
| `--keep-temp` | off | Retain intermediate FASTA/SAM files |
| `--resume` | off | Skip minimap2 if SAM files already exist |

For HPC clusters:

```bash
bash submit_slurm.sh -i <manifest.csv> -r <reference.fa> -a <assembly> -o results/ -t 64
```

---

## Outputs

All outputs are written to `--output-dir`:

| File | Description |
|---|---|
| `{prefix}_remapped_{assembly}.csv` | Full remapped manifest — coordinates + 14 quality columns |
| `matchingSNPs_binary_consistantMapping.{assembly}_map` | **Main output.** Final high-quality marker map (tab-delimited, no header) |
| `{prefix}_remapped_{assembly}.bim` | PLINK BIM format (CHR, SNP, 0, POS, REF, ALT) |
| `_matchingSNPs_binary_consistantMapping.vcf` | Final filtered VCF |
| `_matchingSNPs.vcf` | VCF after design-conflict filter |
| `_matchingSNPs_binary.vcf` | VCF after polymorphic-site filter |
| `allele_usage_decision.txt` | Per-SNP orientation decision (`as_is` / `complement`) |
| `remapping_Report.txt` | Alignment, pair-filtering, position-resolution, and MAPQ summary |
| `QC_Report.txt` | Marker counts at each QC filter stage |
| `remap_assessment/` | MAPQ histograms and known-assembly benchmarks |

### Remapped CSV Quality Columns

The remapped CSV includes 14 new columns beyond the manifest fields:

| Column | Description |
|---|---|
| `Chr_{assembly}` | Chromosome (`"0"` = unmapped/ambiguous) |
| `MapInfo_{assembly}` | **Final 1-based position** (probe-derived or CIGAR-derived; see CoordSource) |
| `Strand_{assembly}` | TopGenomicSeq alignment strand (`+`, `−`, `N/A`) |
| `Ref_{assembly}` / `Alt_{assembly}` | Alleles in alignment strand (strand-normalised by qc_filter.py) |
| `MAPQ_TopGenomicSeq` / `MAPQ_Probe` | Alignment MAPQ scores; MAPQ_Probe=NaN for `topseq_only` markers (no probe alignment); populated normally for `probe_only` |
| `DeltaScore_TopGenomicSeq` | AS gap between best and 2nd-best TopSeq alignments; −1 = uniquely placed |
| `QueryCov_TopGenomicSeq` | Fraction of TopSeq query in aligned (M/=/X) ops |
| `SoftClipFrac_TopGenomicSeq` | Fraction of TopSeq query that is soft-clipped |
| `CoordProbe_{assembly}` | Raw probe-derived coordinate (before any CIGAR override) |
| `Coord_TopSeqCIGAR_{assembly}` | Coordinate from TopSeq CIGAR walk; 0 if SNP target in soft clip |
| `CoordDelta_{assembly}` | `\|CoordProbe − Coord_TopSeqCIGAR\|`; −1 if CIGAR unavailable |
| `CoordSource_{assembly}` | `"probe"` or `"cigar"` — which coordinate is in MapInfo |
| `RefBaseMatch_{assembly}` | Does genome ref base at MapInfo match Ref after strand normalisation? |

### New columns (v2 algorithm)

| Column | Meaning |
|---|---|
| `AlignmentStatus_{assembly}` | Diagnostic census of which alignment sources had hits: `gp1` (both TopSeq alleles + probe), `gp2` (one TopSeq + probe), `gp3` (both TopSeq, no probe), `gp4` (one TopSeq, no probe), `gp5` (probe only), `unmapped`. Computed before any filtering. |
| `anchor_{assembly}` | Which source(s) determined the final coordinate: `topseq_n_probe`, `topseq_only`, `probe_only`, or `N/A` (unmapped/ambiguous). |
| `tie_{assembly}` | How a multi-locus tie was resolved: `unique`, `AS_resolved`, `dAS_resolved`, `NM_resolved`, `CoordDelta_resolved`, `scaffold_resolved`, `ambiguous`, or `N/A`. |
| `RefAltMethodAgreement_{assembly}` | For SNPs: relationship between genome-lookup and NM-based Ref/Alt calls. For indels: whether deletion Ref was confirmed by a genome fetch. See values below. |

**`RefAltMethodAgreement` values for SNPs:**

| Value | Meaning |
|---|---|
| `NM_match` | Genome lookup and NM comparison both succeeded and agree |
| `NM_unmatch` | Both succeeded but disagree — genome result used (flag for QC) |
| `NM_tied` | Genome succeeded; NM was tied — genome result used |
| `NM_N/A` | Genome succeeded; NM not applicable (probe_only marker) |
| `NM_only` | Genome lookup failed; NM result used |
| `ambiguous` | Both methods failed — Chr=0 |

**What is NM?** `NM` is the `NM:i:<n>` edit-distance tag written by minimap2 into each SAM alignment record. It counts mismatches and gap opens between the aligned sequence and the reference — it is **not** derived from CIGAR walking. CIGAR walking in this pipeline is used only for coordinate computation (`parse_cigar_to_ref_pos`, `get_probe_coordinate`). `NM_match`/`NM_unmatch` markers indicate that the genome and NM methods gave the same/different Ref allele assignment; `NM_unmatch` markers are worth inspecting for nearby variants.

### Map File Format

`matchingSNPs_binary_consistantMapping.{assembly}_map` — tab-delimited, no header:

| Column | Description |
|---|---|
| chr | Chromosome |
| pos | Base-pair position |
| snpID | Marker name |
| SNP_alleles | Manifest alleles (e.g. `A,G`) |
| genomic_alleles | + strand alleles matching SNP_alleles order |
| SNP_ref_allele | The SNP allele corresponding to the reference |
| genomic_ref_allele | The reference allele on the + strand |
| allele_usage_decision | `as_is` or `complement` |

---

## QC Filter Cascade

Filters applied sequentially by `qc_filter.py`; `QC_Report.txt` records counts at each stage:

1. **Unmapped** — `Strand == N/A`
2. **MAPQ** — `MAPQ_TopGenomicSeq < --mapq-topseq`
3. **CoordDelta** *(optional)* — `CoordDelta > --coord-delta` or `anchor_{assembly} == "topseq_only"`; disabled by default (`--coord-delta -1`)
4. **Strand agreement** *(optional)* — `StrandAgreementAsExpected == False`; disabled by default (`--require-strand-agreement`)
5. **Design conflict** — SNPs: strand-normalised Ref ≠ genome ref base at MapInfo; deletions: pysam fetch of ref sequence at MapInfo ≠ gref; insertions: always pass
6. **Exclude indels** *(optional)* — remove all indel markers; disabled by default (`--exclude-indels`)
7. **Polymorphic sites** — multiple Ref/Alt assignments at the same Chr:Pos
8. **Consistency** — probe + TopGenomicSeq SAM record count ≠ 3

For Equine80select v2 → EquCab3 (default `--mapq-topseq 30`, no `--coord-delta`):

```
Input markers:                84,319
After unmapped filter:        83,923   (−396)
After MAPQ filter (≥30):      82,406   (−1,517)
After design conflict:        82,178   (−228)
After polymorphic filter:     82,147   (−31)
After consistency filter:     81,491   (−656)
```

With `--coord-delta 0` (removes CoordDelta>0 and all topseq_only):

```
After CoordDelta filter:      81,479   (−927: 186 CoordDelta>0, 741 topseq_only)
Final markers:                81,347
```

---

## Running the Tests

```bash
conda activate remap
pytest tests/ -v
```

Expected: **86 tests** in `tests/test_remap_manifest.py` + **30 tests** in `tests/test_qc_filter.py` + **46 tests** in `tests/test_benchmark_compare.py` (~7 minutes total including integration tests).

The two integration tests in `test_benchmark_compare.py` require real data files in `backup_original/` and `results_E80selv2_to_equCab3/`.

---

## Benchmarking Remapping Accuracy

For assemblies where the true coordinates are known (e.g., EquCab3-native manifest), run `benchmark_compare.py` after each pipeline execution:

```bash
python scripts/benchmark_compare.py \
    --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
    --assembly  equCab3 \
    --output-dir results_E80selv2_to_equCab3/benchmark/
```

To compare probe-derived, CIGAR-derived, and final coordinates in a three-way accuracy table:

```bash
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
    --assembly  equCab3
```

**Current accuracy (v2 manifest → EquCab3, 82,222 benchmarked markers):**
- 98.7% correct (Chr + MapInfo + Strand all match ground truth)
- 99.7% coord-accurate (Chr + MapInfo match, strand may differ)
- 0.2% coord_off · 0.2% unmapped

| Benchmark output | Contents |
|---|---|
| `benchmark_{ts}.tsv` | Per-marker outcome |
| `benchmark_{ts}_report.txt` | Summary + accuracy stratified by CoordDelta |
| `benchmark_{ts}_diff.txt` | Category transitions vs `--baseline` |

---

## Pipeline Methodology

### Dual-Alignment Strategy

For each marker, two sequences are aligned to the reference with `minimap2 -ax sr -N 5`:

1. **TopGenomicSeq** — the full genomic context `PREFIX[AlleleA/AlleleB]SUFFIX` is split into two candidates (one per allele). Both primary and secondary alignments (up to 5) are retained. Triples are ranked by AS score, then ΔAS, then NM. The candidate with lower NM at the winning locus is the reference allele.

2. **Probe** (`AlleleA_ProbeSeq`, 50 bp) — aligned independently. Must map to the same chromosome and overlap the TopSeq alignment window. No strand constraint — bottom-strand probes align opposite to TopGenomicSeq and are still valid.

**Coordinate selection:** For each mapped marker, the probe coordinate and the CIGAR-derived coordinate from TopGenomicSeq are both computed. `CoordDelta = |probe_coord − CIGAR_coord|`. If `CoordDelta ≥ 2`, the CIGAR coordinate is used as `MapInfo` (benchmark shows ~85% accuracy vs ~4% with probe at delta>10). Otherwise probe is used.

**TopSeq-only rescue:** When no valid (TopSeq × probe) triple exists but TopSeq aligned, the SNP coordinate is derived directly from the TopGenomicSeq CIGAR walk. These markers receive `anchor_{assembly} = topseq_only` (72.9% empirical accuracy vs ground truth). In v2→EquCab3, 1,336 markers are rescued this way.

### Infinium Chemistry

- **Infinium II**: variant = base immediately after probe 3′ end
- **Infinium I**: variant = last base of probe
- On minus strand, the probe's physical 3′ end maps to the alignment start position

### Allele Usage Decision

`as_is` or `complement` via XOR of three conditions:
- `IlmnStrand != SourceStrand`
- `SourceSeq != TopGenomicSeq`
- `Strand_{assembly} == '-'`

---

## Optional: Molly Cross-Validation

```bash
bash scripts/compare_molly.sh \
    -b results/Equine80select_remapped_equCab3.bim \
    -m /path/to/MNEc670k.unique_remap.FINAL.csv \
    -o results/molly_comparison/
```

---

## Citation

> Tamer A. Mansour. "A Context-Aware Computational Pipeline for High-Precision Remapping of Genotyping Arrays: Updating the Equine80select Manifest to EquCab3." https://github.com/drtamermansour/Equine80select_remapper, 2025.

## License

MIT License — see [LICENSE](LICENSE).
