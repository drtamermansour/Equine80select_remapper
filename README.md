# Array Manifest Remapper

A pipeline for **moving genotyping array markers from one reference genome assembly to another**.

When a genome assembly is updated, the genomic coordinates of every marker on
an Illumina genotyping array can shift — sometimes by a few base pairs, sometimes
by megabases. Naively lifting coordinates with a chain file misses markers and
loses information about strand, alleles, and confidence. This tool re-aligns
each marker's actual probe sequence and surrounding context to the new
reference and produces a high-confidence map suitable for downstream
genotype-calling, GWAS, and imputation pipelines (PLINK2, Beagle, etc.).

> **Originally built for** the Equine80select array (~82,000 markers) moving
> from EquCab2 to EquCab3 reference assemblies, but works on any Illumina
> Infinium manifest paired with any reference FASTA.

---

## What you get

For every marker in the input manifest, the pipeline produces:

- a **chromosome and base-pair position** on the target reference,
- a **strand** assignment (which DNA strand the alleles live on),
- **Ref / Alt alleles** in VCF-standard orientation,
- a **confidence label** so downstream users can filter low-quality remappings,
- and a per-marker trace showing exactly which quality-control rule (if any)
  removed it from the final output.

The headline output is a tab-separated file —
`{prefix}_allele_map_{assembly}.tsv` — a per-marker crosswalk between the
manifest's original SNP alleles and the new-assembly genomic alleles, plus a
PLINK BIM / final VCF you can feed directly to PLINK2 or Beagle.

---

## Setup

You need [Conda](https://docs.conda.io) (or Mamba) and a reference genome
FASTA file.

```bash
git clone git@github.com:drtamermansour/Equine80select_remapper.git
cd Equine80select_remapper
bash install.sh        # creates a conda environment called 'remap'
conda activate remap
```

The `install.sh` script installs `minimap2`, `samtools`, `bcftools`, `pysam`,
and `pandas` — everything the pipeline needs.

---

## Quick start

```bash
bash run_pipeline.sh \
    -i your_manifest.csv \
    -r reference_genome.fa \
    -a equCab3 \
    -o results/
```

That's it. Three required inputs: the **manifest** (`-i`), the **reference
FASTA** (`-r`), and an **assembly label** (`-a`) used to name the output
columns. Results land in `results/`:

```
results/
├── remapping/                                      ← step 1: realignment
│   ├── {prefix}_remapped_{assembly}.csv            full marker table with quality columns
│   └── remapping_Report.txt                        per-marker decision summary
└── qc/                                             ← step 2: quality filtering
    ├── {prefix}_allele_map_{assembly}.tsv          ★ main output — allele crosswalk
    ├── {prefix}_remapped_{assembly}.bim            PLINK BIM
    ├── {prefix}_remapped_{assembly}.vcf            final filtered VCF
    ├── {prefix}_remapped_{assembly}_traced.csv     per-marker filter trace
    └── QC_Report.txt                               per-stage filter counts
```

For HPC clusters: `bash submit_slurm.sh -i ... -r ... -a ... -o results/ -t 64`.

---

## Common options

The most useful flags. See [docs/cli_reference.md](docs/cli_reference.md) for everything.

| Flag | Default | What it does |
|---|---|---|
| `-t / --threads` | `4` | Threads for minimap2 (use more on HPC) |
| `--preset` | none | One-knob strictness: `strict` / `default` / `permissive`. Tunes the strictness + threshold + include/exclude flags as a bundle; individual flags override. |
| `--min-anchor` | `topseq` | How permissive to be about which markers count: `dual` (strictest — `topseq_n_probe` only), `topseq` (also `topseq_only`), `probe` (most permissive — also `probe_only`) |
| `--include-indels` | off | By default, indel markers are dropped from the final outputs. Pass this to keep them. |
| `--include-ambiguous-snps` | off | By default, A/T and C/G SNPs are dropped (their alleles can't tell strand apart). Pass this to keep them. |
| `--resume` | off | Skip realignment if it has already been done — useful when iterating on filter settings. |

---

## Documentation

Topic-focused guides for users and developers:

| Page | What's in it |
|---|---|
| [docs/algorithm_overview.md](docs/algorithm_overview.md) | How the pipeline decides where each marker goes |
| [docs/decision_tree_simple.md](docs/decision_tree_simple.md) | The full per-marker decision flow as a flowchart |
| [docs/cli_reference.md](docs/cli_reference.md) | Every CLI flag for `run_pipeline.sh`, `remap_manifest.py`, `qc_filter.py` |
| [docs/qc_filters.md](docs/qc_filters.md) | The 11-stage quality-control cascade explained |
| [docs/output_formats.md](docs/output_formats.md) | All output files and their columns |
| [docs/benchmarking.md](docs/benchmarking.md) | How to measure remapping accuracy against a known-good manifest |
| [docs/benchmarking_vs_liftover.md](docs/benchmarking_vs_liftover.md) | Head-to-head benchmark vs. UCSC liftOver and CrossMap |
| [docs/scaffold_filtering.md](docs/scaffold_filtering.md) | Pre-pipeline step to clean unplaced-scaffold haplotypes from the reference |

Deeper dives into specific concepts:

| Page | What's in it |
|---|---|
| [docs/strand_explained.md](docs/strand_explained.md) | The three "strand" columns in an Illumina manifest and what they each mean |
| [docs/CIGAR_walk.md](docs/CIGAR_walk.md) | How a CIGAR string is walked to derive a coordinate |
| [docs/NM_comparison.md](docs/NM_comparison.md) | How alignment edit-distance (`NM`) is used to choose Ref vs Alt for indels |
| [docs/why_we_right.md](docs/why_we_right.md) | Worked examples where the pipeline correctly placed markers that the manifest had wrong |
| [docs/cant_remap.md](docs/cant_remap.md) | Categories of marker the pipeline cannot remap, and why |
| [docs/scaffold_haplotype_thresholds.md](docs/scaffold_haplotype_thresholds.md) | Recommended thresholds for the alt-haplotype scaffold filter |

---

## Running the tests

Unit tests only (no pipeline data needed):

```bash
conda activate remap
pytest tests/ -v -k "not integration"
```

Full suite including the three benchmark-integration tests (require a real
pipeline output directory — typically the `results/` folder from an earlier
run — plus the source Illumina manifest that was fed to the pipeline):

```bash
pytest tests/ -v \
    --results-dir /path/to/results \
    --manifest    /path/to/source_manifest.csv
```

Without these flags, the three integration tests in
`tests/test_benchmark_compare.py` fail fast with a clear
"Integration tests require --results-dir" / "require --manifest" message;
the 292 unit tests run regardless.

---

## Citation

> Tamer A. Mansour. *A Context-Aware Computational Pipeline for High-Precision
> Remapping of Genotyping Arrays: Updating the Equine80select Manifest to
> EquCab3.* https://github.com/drtamermansour/Equine80select_remapper, 2025.
