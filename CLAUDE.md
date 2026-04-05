# CLAUDE.md — Developer Context for Array Manifest Remapper

This file provides context for AI assistants (Claude Code) and human contributors working on this repository. Read it before making any changes.

---

## Project Purpose

This pipeline remaps Illumina genotyping array manifests from one reference genome assembly to another. The primary use case is remapping the **Equine80select** array (~82k SNP markers) from **EquCab2** to **EquCab3**. The core challenge is that simple coordinate lifting fails in repetitive regions and inverted loci, so a dual-alignment approach is used instead.

**The main output used in downstream analysis** is:
```
matchingSNPs_binary_consistantMapping.{assembly}_map
```
This tab-delimited file maps ~80,197 high-quality markers and is consumed directly by GWAS/imputation pipelines (e.g., PLINK2, Beagle).

---

## Architecture

```
install.sh
    └── Creates conda env 'remap' (minimap2, samtools, bcftools, pysam, pandas)

run_pipeline.sh                          ← main entry point (bash, arg-parsed)
    ├── samtools faidx + vcf_contigs.txt ← reference prep (step 1)
    ├── scripts/remap_manifest.py        ← dual alignment + coordinate resolution (step 2)
    └── scripts/qc_filter.py            ← QC cascade + all output generation (step 3)

scripts/compare_molly.sh                 ← optional standalone cross-validation
submit_slurm.sh                          ← SLURM wrapper for HPC
```

---

## Module Details

### `scripts/remap_manifest.py`

**Input:** Illumina manifest CSV, reference FASTA
**Output:** Remapped manifest CSV with 7 new columns (see below)

Key steps:
1. Parses the `[Assay]` section of the Illumina manifest (skipping header/footer)
2. Builds two FASTA files: `temp_probes.fasta` (50bp AlleleA probes) and `temp_topseq.fasta` (two sequences per marker: `PREFIX+AlleleA+SUFFIX` and `PREFIX+AlleleB+SUFFIX`)
3. Runs `minimap2 -ax sr` on both — `-N 5` for probes to allow secondary alignments for overlap checking
4. Parses TopGenomicSeq SAM: lower NM (edit distance) candidate = reference allele
5. Selects best probe alignment by overlap with TopGenomicSeq window, MAPQ as tiebreaker
6. Calculates variant position via `get_probe_coordinate()` (probe-based, primary) or `parse_cigar_to_ref_pos()` (TopSeq CIGAR fallback)

**CLI flags:** `-i`, `-r`, `-o`, `-a` (assembly name), `--threads`, `--temp-dir`, `--keep-temp`, `--resume`

### `scripts/qc_filter.py`

**Input:** Remapped manifest CSV + reference FASTA + vcf_contigs.txt + SAM files (for consistency check)
**Output:** All filtered VCFs, BIM, final map file, QC_Report.txt, remap_assessment/

Filter stages (in order):
1. `Strand_{assembly} == 'N/A'` → unmapped, remove
2. `MAPQ_TopGenomicSeq < --mapq-topseq` → low-confidence mapping, remove
3. Remapped Ref ≠ genome Ref at that position → design conflict, remove
4. Multiple Ref/Alt combinations at same Chr:Pos → polymorphic site, remove
5. Count of (topseq_A + topseq_B + probe) SAM records ≠ 3 → inconsistent, remove

Also computes:
- `allele_usage_decision.txt` (as_is / complement per marker)
- MAPQ histograms in `remap_assessment/`
- Known-assembly benchmark (markers where original assembly == target, checks coordinate match)

**CLI flags:** `-i`, `-r`, `-v`, `-a`, `-o`, `--mapq-topseq`, `--mapq-probe`, `--temp-dir`, `--prefix`, `--topseq-sam`, `--probe-sam`

---

## Key Output Column Names

The assembly name (from `--assembly` / `-a`) is embedded in column names. For EquCab3:

| Column | Type | Meaning |
|---|---|---|
| `Chr_equCab3` | str | Chromosome on new assembly (`"0"` = unmapped) |
| `MapInfo_equCab3` | int | 1-based base-pair position |
| `Strand_equCab3` | str | `+`, `-`, or `N/A` (from TopGenomicSeq alignment) |
| `Ref_equCab3` | str | Reference allele (from TopGenomicSeq, NOT strand-normalized) |
| `Alt_equCab3` | str | Alternate allele (same) |
| `MAPQ_TopGenomicSeq` | int | MAPQ of winning TopGenomicSeq alignment |
| `MAPQ_Probe` | int | MAPQ of selected probe alignment; **0 = TopSeq fallback used** |

**Important:** `Ref_equCab3` and `Alt_equCab3` are in the **alignment strand** of the TopGenomicSeq winner, not necessarily the forward (+) strand. The `qc_filter.py` strand-normalizes them before VCF/map output using `strand_normalize(allele, strand)`.

---

## Key Biological Concepts

### Illumina Manifest Format
- **`AlleleA_ProbeSeq`**: the 50bp physical probe sequence (Infinium I: two probes; Infinium II: one probe)
- **`AlleleB_ProbeSeq`**: present for Infinium I (two different probes), absent/NaN for Infinium II
- **`TopGenomicSeq`**: longer genomic context `PREFIX[AlleleA/AlleleB]SUFFIX` — the square-bracket section is the SNP
- **`SNP`**: the marker's allele pair in Illumina notation, e.g. `[A/G]`
- **`IlmnStrand`**, **`SourceStrand`**: describe the relationship between the probe and the reference strand (TOP/BOT or PLUS/MINUS)

### Infinium Chemistry
- **Infinium II**: single probe, extends into the SNP site — variant position is 1 base AFTER the probe's 3' end
- **Infinium I**: two probes (one per allele), both ending AT the SNP — variant position IS the probe's last base
- Determined in code by: `pd.isna(row['AlleleB_ProbeSeq'])` → True = Infinium II

### Dual-Alignment Strategy
TopGenomicSeq is the **strand authority** — it determines chromosome, strand, and Ref/Alt. The probe provides **coordinate precision** when it aligns consistently. If probe alignment disagrees with TopGenomicSeq strand or falls on a different chromosome, the coordinate falls back to parsing the TopGenomicSeq CIGAR string.

### Allele Usage Decision (8-condition XOR logic)
Controls how manifest SNP alleles (`[A/B]` column) map to forward-strand genomic alleles:

```
flip = (IlmnStrand != SourceStrand) XOR (SourceSeq != TopGenomicSeq) XOR (Strand == '-')
decision = 'complement' if flip else 'as_is'
```

Indel markers get the prefix `indel_` (e.g., `indel_as_is`). These are tracked but excluded from VCF output.

---

## Final Map File Format

`matchingSNPs_binary_consistantMapping.{assembly}_map` — tab-delimited, no header, 8 columns:

```
chr  pos  snpID  SNP_alleles  genomic_alleles  SNP_ref_allele  genomic_ref_allele  decision
```

- `SNP_alleles`: allele pair from manifest `[A/B]` column, comma-separated (e.g., `A,G`)
- `genomic_alleles`: forward-strand alleles in the **same order** as SNP_alleles (e.g., `G,A` means SNP allele A ↔ genomic G)
- `SNP_ref_allele`: the manifest allele corresponding to the genomic reference
- `genomic_ref_allele`: the reference allele on the forward strand
- `decision`: `as_is` or `complement` (or `indel_*` variant)

---

## Common Pitfalls for Contributors

1. **Column names depend on `--assembly` flag.** If you pass `-a equCab3`, columns are `Chr_equCab3`. If you pass `-a EquCab3`, they are `Chr_EquCab3`. These must match between `remap_manifest.py` and `qc_filter.py` — `run_pipeline.sh` handles this automatically.

2. **Resume after step 3 failure.** If `qc_filter.py` fails after `remap_manifest.py` has already completed, re-run with `--resume` to skip the expensive minimap2 alignment step. Temp files (including SAM files) are never deleted until the full pipeline completes successfully, so `--resume` always has access to the SAM files needed for the consistency filter — no special flags required on the original run.

3. **Consistency check requires SAM files.** `qc_filter.py` needs `temp_topseq.sam` and `temp_probe.sam` for the consistency filter (step 5). When run via `run_pipeline.sh`, these files always exist at step 3 (they are only deleted at the very end of a successful run). If calling `qc_filter.py` directly on the output of an already-completed pipeline (where temp files were already cleaned up), pass `--topseq-sam` and `--probe-sam` explicitly, or the consistency filter is skipped with a warning.

4. **bcftools must be on PATH.** `qc_filter.py` calls `bcftools norm` as a subprocess. Activating the `remap` conda environment provides this.

5. **minimap2 `-ax sr` mode** is designed for short genomic reads (~50–300 bp). Both probes (~50 bp) and TopGenomicSeq (~150 bp) fall in this range.

6. **The `Chr` column in the remapped CSV is a string**, not an integer — chromosomes like `Un_NW_*` are valid. Always use `dtype={col_chr: str}` when loading.

7. **Deletion correction on minus strand.** When a marker is a deletion (len(Ref) > len(Alt)) and the strand is `-`, `remap_manifest.py` subtracts `del_len` from the raw probe coordinate. Do not remove this logic — it reflects the directional asymmetry of probe 3' end positioning.

---

## How to Run (Equine80select → EquCab3)

```bash
conda activate remap

bash run_pipeline.sh \
    -i backup_original/Equine80select_24_20067593_B1.csv \
    -r equCab3/equCab3_genome.fa \
    -a equCab3 \
    -o results/ \
    -t 8
```

Expected: ~80,197 markers in `results/matchingSNPs_binary_consistantMapping.equCab3_map`

---

## How to Add a New Assembly

No code changes are needed. Just:

1. Download the target assembly FASTA
2. Strip `chr` prefix if present: `sed 's/>chr/>/' raw.fa > assembly.fa`
3. Index it: `samtools faidx assembly.fa`
4. Run the pipeline with the new `-r` and `-a` values:
   ```bash
   bash run_pipeline.sh -i manifest.csv -r /path/to/assembly.fa -a assemblyName -o results/
   ```

Output column names and file names will use your `assemblyName` automatically.

---

## File Inventory

### Scripts (production)
| File | Purpose |
|---|---|
| `install.sh` | One-time conda environment setup |
| `run_pipeline.sh` | Main orchestrator (arg-parsed, calls both Python scripts) |
| `submit_slurm.sh` | SLURM wrapper for HPC submission |
| `myRun.sh` | Project-specific SLURM run script for Equine80select → EquCab3; records the exact flags used (`--assembly equCab3 --threads 64 --keep-temp --resume --mapq-topseq 1`) |
| `scripts/remap_manifest.py` | Core remapping: dual alignment + coordinate resolution |
| `scripts/qc_filter.py` | QC cascade, allele decision, VCF/BIM/map generation |
| `scripts/scaffold_haplotype_analyzer.py` | Pre-pipeline: align all unplaced scaffolds to chromosomal reference. |
| `scripts/filter_scaffold_haplotypes.py` | Pre-pipeline: filter the aligned unplaced scaffolds, output `scaffold_haplotype_analysis/alt_haplotype_candidates.tsv` |
| `scripts/compare_molly.sh` | Optional: cross-validation vs. MNEc670k dataset |

### Reference data
| File/Dir | Purpose |
|---|---|
| `equCab3/equCab3_genome.fa` | EquCab3 reference genome FASTA |
| `equCab3/equCab3_genome.fa.fai` | samtools index |
| `equCab3/vcf_contigs.txt` | VCF header contig definitions (generated by pipeline) |
| `equCab2/equCab2_genome.fa` | EquCab2 reference (symlink; needed for EquCab2 benchmarking) |
| `backup_original/` | Original Illumina manifest CSV (do not modify) |

### Generated outputs of remapping (from completed Equine80selectv1 → EquCab3 run)
| File | Description |
|---|---|
| `Equine80select_24_20067593_B1_remapped.csv` | Remapped manifest (all markers) |
| `matchingSNPs_binary_consistantMapping.EquCab3_map` | **Main output** (80,197 markers) |
| `Equine80select_remapped_equCab3.bim` | PLINK BIM (80,223 markers) |
| `_matchingSNPs_binary_consistantMapping.vcf` | Final VCF |
| `allele_usage_decision.txt` | Allele orientation decisions |
| `MNEc670k_remap.tab` | Molly dataset (cleaned tab form) |
| `molly_remapped_*.txt` | Cross-validation results |
| `remap_assessment/` | MAPQ histograms, benchmark mismatches |
| `inconsistent_remapped.csv` | Markers excluded by consistency filter |
| `polymorphic_positions.txt` | Positions excluded by polymorphic filter |
| `temp_probes.fasta`, `temp_topseq.fasta` | Alignment input sequences (kept for debugging) |
| `temp_probe.sam`, `temp_topseq.sam` | Alignment outputs (needed for consistency check) |
| `_tmp.*` | Intermediate allele tables (generated by pipeline) |

### Scratch / archive (ignore)
| File/Dir | Purpose |
|---|---|
| `temp/` | Archive folder — not part of the pipeline. Contains `master_pipeline.sh` (superseded by `run_pipeline.sh`) and two troubleshooting script fragments (`tmp_benchmark.sh`, `tmp_test_remapping.sh`). **Ignore for development.** |
