# Array Manifest Remapper

A computational pipeline for remapping Illumina genotyping array manifests between reference genome assemblies. It uses a **context-aware dual-alignment strategy** — aligning both the short physical probe (50 bp) and the longer `TopGenomicSeq` context sequence — to ensure high-fidelity coordinate conversion, correct strand assignment, and precise Ref/Alt allele determination consistent with VCF standards.

Originally developed to remap the **Equine80select** array from EquCab2 to EquCab3.

## Prerequisites

- **Conda** or **Mamba** package manager
- **Reference genome FASTA** for the target assembly (indexed with `samtools faidx`)

## Setup

```bash
git clone https://github.com/drtamermansour/Equine80select_remapper.git
cd Equine80select_remapper

# Create the 'remap' conda environment with all dependencies
bash install.sh
conda activate remap
```

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
| `--keep-temp` | off | Retain intermediate FASTA/SAM files |
| `--resume` | off | Skip step 2 if the remapped CSV already exists (useful after a step 3 failure) |

For HPC clusters, use the SLURM wrapper:

```bash
bash submit_slurm.sh -i <manifest.csv> -r <reference.fa> -a <assembly> -o results/ -t 64
```

## Outputs

All outputs are written to `--output-dir`:

| File | Description |
|---|---|
| `{prefix}_remapped_{assembly}.csv` | Full remapped manifest with new coordinates and MAPQ scores |
| `matchingSNPs_binary_consistantMapping.{assembly}_map` | **Main output.** Final high-quality marker map |
| `{prefix}_remapped_{assembly}.bim` | PLINK BIM format (CHR, SNP, 0, POS, REF, ALT) |
| `_matchingSNPs_binary_consistantMapping.vcf` | Final filtered VCF |
| `_matchingSNPs.vcf` | VCF after design-conflict filter |
| `_matchingSNPs_binary.vcf` | VCF after polymorphic-site filter |
| `allele_usage_decision.txt` | Per-SNP orientation decision (as_is / complement) |
| `QC_Report.txt` | Marker counts at each filter stage |
| `remap_assessment/` | MAPQ histograms and known-assembly benchmarks |

### Map File Format

`matchingSNPs_binary_consistantMapping.{assembly}_map` — tab-delimited, no header:

| Column | Description |
|---|---|
| chr | Chromosome on the new assembly |
| pos | Base-pair position |
| snpID | Marker name |
| SNP_alleles | Manifest alleles (e.g. `A,G`) |
| genomic_alleles | + strand alleles matching SNP_alleles order |
| SNP_ref_allele | The SNP allele corresponding to the reference |
| genomic_ref_allele | The reference allele on the + strand |
| allele_usage_decision | `as_is` or `complement` |

## QC Filter Cascade

Markers are removed progressively; `QC_Report.txt` records counts at each stage:

1. **Unmapped** — `Strand == N/A` (failed TopGenomicSeq alignment)
2. **MAPQ** — `MAPQ_TopGenomicSeq < --mapq-topseq`
3. **Design conflict** — remapped Ref allele does not match genome Ref at that position
4. **Polymorphic sites** — position shared by multiple markers with different Ref/Alt assignments
5. **Consistency** — probe + TopGenomicSeq alignment count ≠ 3

For Equine80select → EquCab3 (default settings):
```
Input markers:                81,974
After unmapped filter:        81,945   (-29)
After MAPQ filter (>=30):     ~81,9xx
After design conflict:        81,672   (-273)
After polymorphic filter:     81,660   (-12)
After consistency filter:     80,197   (-1,463)
```

## Optional: Molly Cross-Validation

For Equine80select, results can be cross-validated against the MNEc670k Molly remapping:

```bash
bash scripts/compare_molly.sh \
    -b results/Equine80select_remapped_equCab3.bim \
    -m /path/to/MNEc670k.unique_remap.FINAL.csv \
    -o results/molly_comparison/
```

## Pipeline Methodology

### Dual-Alignment Strategy

For each marker, two sequences are aligned to the reference with `minimap2 -ax sr`:

1. **TopGenomicSeq** — the full genomic context `PREFIX[AlleleA/AlleleB]SUFFIX` is split into two candidates (one per allele). The candidate with lower edit distance (NM tag) is the reference allele. This alignment determines chromosome, strand, and Ref/Alt.

2. **Probe** (`AlleleA_ProbeSeq`, 50 bp) — if the probe aligns to the same chromosome and overlaps the TopGenomicSeq window on the same strand, its 3' end position is used for high-precision coordinate calculation. Otherwise, the coordinate falls back to parsing the TopGenomicSeq CIGAR string.

### Infinium Chemistry Handling

- **Infinium II**: variant is the base immediately after the probe 3' end
- **Infinium I**: variant is the last base of the probe
- On the minus strand, the probe's physical 3' end maps to the alignment start position

### Allele Usage Decision

Each marker gets an `as_is` or `complement` decision, controlling how the manifest's SNP alleles (`[A/B]` column) map to forward-strand genomic alleles. The decision is determined by three binary conditions (XOR logic):
- `IlmnStrand != SourceStrand`
- `SourceSeq != TopGenomicSeq`
- `Strand_{assembly} == '-'`

## Citation

If you use this pipeline in your research, please cite:

> Tamer A. Mansour. "A Context-Aware Computational Pipeline for High-Precision Remapping of Genotyping Arrays: Updating the Equine80select Manifest to EquCab3." https://github.com/drtamermansour/Equine80select_remapper, 2025.

## License

MIT License — see [LICENSE](LICENSE).
