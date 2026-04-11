# CLAUDE.md — Developer Context for Array Manifest Remapper

Read this before making any changes. It is the single source of truth for how the pipeline works.

---

## Project Purpose

Remaps Illumina genotyping array manifests between reference genome assemblies using dual alignment (probe + TopGenomicSeq via minimap2). Primary use case: **Equine80select** array (~82k SNPs) remapped from EquCab2 → EquCab3.

**Main output used downstream:**
```
matchingSNPs_binary_consistantMapping.{assembly}_map
```
Consumed directly by GWAS/imputation pipelines (PLINK2, Beagle).

---

## Architecture

```
install.sh           → conda env 'remap' (minimap2, samtools, bcftools, pysam, pandas)
run_pipeline.sh      → main orchestrator; calls steps 1–3
  step 1: samtools faidx + vcf_contigs.txt
  step 2: scripts/remap_manifest.py
  step 3: scripts/qc_filter.py
submit_slurm.sh      → SLURM wrapper
myRun.sh             → exact flags used for Equine80select → EquCab3
scripts/benchmark_compare.py  → standalone accuracy benchmark (run manually)
scripts/benchmark_cigar_vs_probe.py  → three-way coord comparison (run manually)
tests/               → pytest suite (106 tests)
```

---

## Pre-pipeline: Scaffold Haplotype Exclusion (Optional)

Modern assemblies include unplaced scaffolds that are often alt-haplotypes of placed chromosomes. Including them causes multi-mapping and reduces confidence. Three scripts handle this **before** the main pipeline:

```bash
# 1. Characterise scaffolds (minimap2 asm5 against placed chromosomes)
python scripts/scaffold_haplotype_analyzer.py -r equCab3/equCab3_genome.fa \
    -o remap_assessment/scaffold_haplotype_analysis/ --threads 8
# → scaffold_summary.tsv with identity_pct, query_coverage_pct, span_to_scaffold_ratio, max_mapq, n_alignment_blocks

# 2. Filter to alt-haplotype candidates (default: Tier 1 thresholds)
python scripts/filter_scaffold_haplotypes.py \
    -i remap_assessment/scaffold_haplotype_analysis/scaffold_summary.tsv \
    -o remap_assessment/scaffold_haplotype_analysis/alt_haplotype_candidates.tsv

# 3. Build cleaned reference (removes identified scaffolds, indexes output)
python scripts/exclude_alt_haplotypes.py \
    --scaffolds remap_assessment/scaffold_haplotype_analysis/alt_haplotype_candidates.tsv \
    --reference equCab3/equCab3_genome.fa --output-dir equCab3_cleaned/
# → {stem}_no_alt_haplotypes.fa, .fai, exclusion_report.txt
```

See `docs/scaffold_haplotype_thresholds.md` for threshold tiers and rationale.

---

## `scripts/remap_manifest.py`

**Input:** Illumina manifest CSV, reference FASTA
**Output:** `{prefix}_remapped_{assembly}.csv` — all input rows plus 17 new columns

### Algorithm

1. Parse `[Assay]` section of manifest.
2. Build `temp_probes.fasta` (50 bp AlleleA probes) and `temp_topseq.fasta` (two sequences per marker: PREFIX+AlleleA+SUFFIX and PREFIX+AlleleB+SUFFIX).
3. Run `minimap2 -ax sr -N 5` on both. **Both primary and secondary alignments are retained.**
4. For each marker, enumerate all valid (TopSeq_allele × TopSeq_align × probe_align) triples via `build_valid_triples`. A triple is valid if:
   - TopSeq chr == probe chr
   - `compute_probe_strand_agreement` returns `"True"` (strand check is a hard filter; `"N/A"` is treated as valid)
   - Among strand-valid probes on the same chromosome, only the highest-overlap one is kept per ts_align.
   - overlap > 0.
5. Rank triples by AS_sum (ts.AS + pb.AS) → ΔAS_sum → NM_sum → CoordDelta → scaffold_resolved → ambiguous. MAPQ is no longer used for ranking (reported as diagnostic only). If no valid triples exist, fall through to `best_topseq_rescue` (anchor=topseq_only) then `best_probe_rescue` (anchor=probe_only).
6. Determine Ref/Alt via `determine_ref_alt_v2`: genome lookup (`resolve_ref_from_genome`, now strand-aware) is the primary method for SNPs; NM comparison runs in parallel. For indels, NM comparison determines Ref/Alt; deletion Ref is validated against the genome (replaces `check_deletion_ref_match` in qc_filter.py). Agreement reported in `RefAltMethodAgreement_{assembly}`.
7. Calculate variant position via `get_probe_coordinate()` using the **probe's own alignment strand** (NOT the TopSeq strand — mixing these caused a ±51 bp bug that is now fixed).
8. **CIGAR cross-validation:** compute `Coord_TopSeqCIGAR` by walking the TopSeq CIGAR to the SNP position. `CoordDelta = |probe_coord − CIGAR_coord|`. If `CoordDelta ≥ 2`, use CIGAR coord as `MapInfo` (`CoordSource = "cigar"`); otherwise use probe coord (`CoordSource = "probe"`). Benchmark accuracy: probe 98.0%, CIGAR 98.6%, final 98.7%.

**topseq_only rescue:** When no valid triple exists but TopSeq aligned, derive SNP coordinate from CIGAR walk on the best TopSeq alignment (highest MAPQ/AS across both alleles). Ref/Alt determined by NM comparison. 1,336/1,394 no-valid-triple markers rescued (95.8%). Failures: 42 markers with SNP in soft-clipped region (CIGAR coord unavailable), 16 with NM tie. `MAPQ_Probe = NaN` and `CoordDelta = −1` for all topseq_only markers.

### Output columns (assembly = `equCab3` example)

| Column | Meaning |
|---|---|
| `Chr_equCab3` | Chromosome (`"0"` = unmapped or ambiguous) |
| `MapInfo_equCab3` | Final 1-based position (probe-derived or CIGAR-derived; see CoordSource) |
| `Strand_equCab3` | `+`, `-`, or `N/A` (TopGenomicSeq alignment strand) |
| `Ref_equCab3` / `Alt_equCab3` | Alleles in alignment strand (NOT strand-normalised; qc_filter.py normalises) |
| `MAPQ_TopGenomicSeq` / `MAPQ_Probe` | Alignment MAPQ scores; MAPQ_Probe=NaN for topseq_only markers (no probe alignment); populated normally for probe_only |
| `AlignmentStatus_equCab3` | Raw alignment census: gp1–gp5 or unmapped (computed before filtering) |
| `anchor_equCab3` | Which source was used for the coordinate: topseq_n_probe / topseq_only / probe_only / N/A |
| `tie_equCab3` | How the winning locus was selected: unique / AS_resolved / dAS_resolved / NM_resolved / CoordDelta_resolved / scaffold_resolved / ambiguous / N/A |
| `RefAltMethodAgreement_equCab3` | Agreement between genome-lookup and NM-based Ref/Alt (see docs/pipeline_decision_tree.md §4 for value definitions) |
| `DeltaScore_TopGenomicSeq` | AS gap between best and 2nd-best TopSeq alignments; −1 if fewer than 2 |
| `QueryCov_TopGenomicSeq` | Fraction of TopSeq query in M/=/X aligned ops; 0.0 for unmapped |
| `SoftClipFrac_TopGenomicSeq` | Fraction of TopSeq query bases soft-clipped; 0.0 for unmapped |
| `CoordProbe_equCab3` | Raw probe-derived coordinate before any CIGAR override |
| `Coord_TopSeqCIGAR_equCab3` | Coordinate from TopSeq CIGAR walk; 0 if unmapped or SNP in soft clip |
| `CoordDelta_equCab3` | `\|CoordProbe − Coord_TopSeqCIGAR\|`; −1 if CIGAR unavailable |
| `CoordSource_equCab3` | `"probe"` or `"cigar"` — which coordinate is in MapInfo |
| `RefBaseMatch_equCab3` | `"True"` / `"False"` / `"N/A"` — genome ref base at MapInfo matches Ref (strand-normalised)? |

**Note:** The v1 algorithm counts (mapped: 82,482 | ref_resolved: 105 | topseq_only: 1,336 | ambiguous: 146 | unmapped: 250) are stale and will differ under the current algorithm. See `remapping/remapping_Report.txt` in the output directory for up-to-date counts.

### Supplementary outputs (in `remapping/` subfolder alongside main CSV)

- `ambiguous_markers.csv`, `scaffold_resolved_markers.csv`, `nm_position_resolved_markers.csv`
- `remapping_Report.txt` — decision summary

### CLI flags
`-i` (manifest), `-r` (reference FASTA), `-o` (output dir), `-a` (assembly name), `--threads`, `--temp-dir`, `--resume`

---

## `scripts/qc_filter.py`

**Input:** Remapped CSV + reference FASTA + `vcf_contigs.txt` + SAM files
**Output:** Filtered VCFs, BIM, final map, `QC_Report.txt`, `remap_assessment/` — all written to the `qc/` subfolder of `--output-dir`

### Filter cascade

1. `Strand_{assembly} == 'N/A'` → unmapped, remove
2. `MAPQ_TopGenomicSeq < --mapq-topseq` → low-confidence, remove
3. *(optional)* `CoordDelta > --coord-delta` OR `anchor_{assembly} == "topseq_only"` → remove; disabled by default (`--coord-delta -1`)
4. *(optional)* `StrandAgreementAsExpected == False` → remove; disabled by default (`--require-strand-agreement`)
5. Design conflict:
   - **SNPs**: strand-normalised Ref ≠ single-base genome Ref from bcftools → remove
   - **Indels (deletion ref)**: pysam fetch of `len(gref)` bases at MapInfo ≠ gref → remove
   - **Indels (insertion ref, `gref == ""`)**: always pass (nothing to verify against reference)
6. *(optional)* `--exclude-indels` → remove all indel markers from output; disabled by default
7. Multiple Ref/Alt assignments at same Chr:Pos → polymorphic, remove
8. Count of SAM records (topseq_A + topseq_B + probe) ≠ 3 → inconsistent, remove

**VCF and BIM indel encoding:** `pos = MapInfo − 1`; `REF = anchor + gref`; `ALT = anchor + galt` (anchor base fetched from FASTA at `MapInfo − 1`; deletion has `galt = ""`; insertion has `gref = ""`). Consistency filter still applies to indels — count is 3 (topseq_A + topseq_B + probe) because all indel markers are Infinium II with one probe.

### Expected counts (v2 manifest → EquCab3, default `--mapq-topseq 30`, no `--coord-delta`)
```
Input:                    84,319
After unmapped:           83,923   (-396)
After MAPQ (>=30):        82,406   (-1,517)
After design conflict:    82,178   (-228)
After polymorphic:        82,147   (-31)
After consistency:        81,491   (-656)
```

With `--coord-delta 0` (removes CoordDelta>0 and all topseq_only):
```
After CoordDelta filter:  81,479   (-927: 186 CoordDelta>0, 741 topseq_only)
Final markers:            81,347
```

### CLI flags
`-i`, `-r`, `-v` (vcf_contigs), `-a`, `-o`, `--mapq-topseq`, `--mapq-probe`, `--coord-delta` (default -1), `--require-strand-agreement`, `--exclude-indels` (default: indels included), `--temp-dir`, `--prefix`, `--topseq-sam`, `--probe-sam`

---

## `scripts/benchmark_compare.py`

Standalone accuracy benchmark for EquCab3-native manifests. Compares remapped Chr/MapInfo/Strand against ground-truth manifest values. Outputs per-marker TSV, summary report with accuracy stratified by CoordDelta, and optional diff vs `--baseline`.

```bash
python scripts/benchmark_compare.py \
    --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
    --assembly  equCab3 --output-dir results_E80selv2_to_equCab3/qc/benchmark/
```

Strand ground truth: `SourceStrand` column (TOP/PLUS→+, BOT/MINUS→-). **Do not use `RefStrand`** — it encodes probe design convention, not alignment strand.

## `scripts/benchmark_cigar_vs_probe.py`

Three-way comparison of `CoordProbe_{assembly}`, `Coord_TopSeqCIGAR_{assembly}`, and `MapInfo_{assembly}` against ground truth. Shows per-CoordDelta-bucket accuracy for each source. Current accuracy: probe 98.0%, CIGAR 98.6%, final 98.7% (82,222 benchmarked markers).

---

## Tests

```bash
conda activate remap
pytest tests/ -v        # 159 tests (~7 min including integration tests)
```

- `tests/test_remap_manifest.py` — 86 tests: `parse_topseq_sam`, `is_placed_chromosome`, `select_best_pair`, `determine_ref_alt`, `parse_cigar_to_ref_pos`, `compute_qcov`, `compute_soft_clip_frac`, `_get_as`, `reverse_complement`, `probe_topseq_orientation`, `compute_probe_strand_agreement`, `extract_candidates`, CIGAR-override-for-indels contract
- `tests/test_qc_filter.py` — 30 tests: `apply_probe_mapq_filter`, `apply_strand_agreement_filter`, `strand_normalize`, `check_deletion_ref_match`, `make_anchor_alleles`, `apply_exclude_indels_filter`
- `tests/test_benchmark_compare.py` — 46 tests: unit + integration for `benchmark_compare.py`; two integration tests require real data in `backup_original/` and `results_E80selv2_to_equCab3/`

---

## Key Biological Concepts

### Illumina manifest columns used by the pipeline

| Column | Used for |
|---|---|
| `AlleleA_ProbeSeq` | 50 bp probe; `AlleleB_ProbeSeq` is NaN for Infinium II |
| `TopGenomicSeq` | Genomic context `PREFIX[A/B]SUFFIX`; split into two alignment candidates |
| `IlmnStrand` | TOP/BOT — used in probe strand agreement check (`remap_manifest.py`) |
| `SourceStrand` | TOP/BOT/PLUS/MINUS — used in probe strand agreement check AND as strand ground truth in benchmark |
| `RefStrand` | +/- — NOT used by the pipeline; encodes probe design convention, not alignment strand |

### Infinium chemistry
- **Infinium II**: `AlleleB_ProbeSeq` is NaN. Variant = base immediately AFTER probe 3' end.
- **Infinium I**: two probes, both ending AT the SNP. Variant = last base of probe.
- On the minus strand, the probe's physical 3' end is at the alignment *start* position.

---

## Final Map File Format

`matchingSNPs_binary_consistantMapping.{assembly}_map` — tab-delimited, no header, 8 columns:
```
chr  pos  snpID  SNP_alleles  genomic_alleles  SNP_ref_allele  genomic_ref_allele  decision
```
- `genomic_alleles`: forward-strand alleles in the **same order** as `SNP_alleles`
- `decision`: `as_is` or `complement` (or `indel_*`) — inferred by matching SNP alleles to remapped genomic alleles

---

## Critical Pitfalls

1. **Column names embed the assembly name.** `-a equCab3` → `Chr_equCab3`. Must match between `remap_manifest.py` and `qc_filter.py`. `run_pipeline.sh` handles this automatically.

2. **Resume after step 3 failure.** Re-run `run_pipeline.sh` with `--resume` to skip minimap2. Temp SAM files are preserved until the full pipeline completes successfully.

3. **Consistency filter needs SAM files.** If calling `qc_filter.py` directly after a completed run (temp files cleaned up), pass `--topseq-sam` and `--probe-sam` explicitly or the filter is skipped with a warning.

4. **Chr column is always a string.** Use `dtype={col_chr: str}` when loading CSVs. Unplaced contigs like `Un_NW_*` are valid chromosome values.

5. **Deletion correction on minus strand.** When Ref is longer than Alt and strand is `-`, `remap_manifest.py` subtracts `del_len` from the raw probe coordinate. This is intentional — do not remove it.

6. **bcftools must be on PATH.** `qc_filter.py` calls `bcftools norm` as a subprocess. `conda activate remap` provides it.

7. **Use `--resume` during development.** Existing SAM files live in `results_E80selv2_to_equCab3/temp/`. Pass `--resume --temp-dir results_E80selv2_to_equCab3/temp` to skip the multi-hour minimap2 alignment step when iterating on parsing or coordinate logic.

8. **CoordDelta on minus-strand indels may be inflated** by `allele_len − 1`. This is expected behaviour; the CIGAR and probe coordinates shift differently for deletions on the minus strand.

9. **`resolve_ref_from_genome` requires `strand`.** The function now complements allele characters before comparing against the forward-strand genome base when `strand == "-"`. All call sites must pass the alignment strand. Omitting it caused silent `None` returns for minus-strand NM ties.

10. **CoordDelta filter uses `anchor_{assembly}`, not `MappingStatus_{assembly}`.** `qc_filter.py --coord-delta` now reads `anchor_{assembly} == "topseq_only"` to identify markers to exclude. Old CSVs without `anchor_{assembly}` will skip the filter with a warning.

---

## File Inventory

### Production scripts
| File | Purpose |
|---|---|
| `install.sh` | One-time conda environment setup |
| `run_pipeline.sh` | Main orchestrator |
| `submit_slurm.sh` | SLURM wrapper |
| `myRun.sh` | Exact flags used for Equine80select → EquCab3 |
| `scripts/remap_manifest.py` | Dual alignment + coordinate resolution → remapped CSV (17 new columns) |
| `scripts/qc_filter.py` | QC cascade, allele decision, VCF/BIM/map generation |
| `scripts/benchmark_compare.py` | Standalone accuracy benchmark (run manually after pipeline) |
| `scripts/benchmark_cigar_vs_probe.py` | Three-way coord comparison: probe vs CIGAR vs final (run manually) |
| `scripts/scaffold_haplotype_analyzer.py` | Pre-pipeline: characterise unplaced scaffolds vs placed chromosomes |
| `scripts/filter_scaffold_haplotypes.py` | Pre-pipeline: filter to alt-haplotype candidates |
| `scripts/exclude_alt_haplotypes.py` | Pre-pipeline: build cleaned reference excluding alt-haplotype scaffolds |
| `scripts/compare_molly.sh` | Optional: cross-validation vs. MNEc670k dataset |

### Reference data
| File/Dir | Purpose |
|---|---|
| `equCab3/equCab3_genome.fa` | EquCab3 reference genome FASTA |
| `equCab2/equCab2_genome.fa` | EquCab2 reference (symlink) |
| `backup_original/` | Original Illumina manifest CSVs (do not modify) |

### Ignore
| File/Dir | Note |
|---|---|
| `temp/` | Archive — superseded scripts; ignore for development |
