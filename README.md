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
| `--require-strand-agreement` | off | Remove markers where probe strand disagrees with expected orientation |
| `--keep-temp` | off | Retain intermediate FASTA/SAM files |
| `--resume` | off | Skip minimap2 if SAM files already exist |

For HPC clusters:

```bash
bash submit_slurm.sh -i <manifest.csv> -r <reference.fa> -a <assembly> -o results/ -t 64
```

---

## Outputs

Outputs are organised into two subdirectories inside `--output-dir`:

```
output-dir/
├── temp/                                        ← intermediate FASTA/SAM files (removed after pipeline)
├── remapping/                                   ← remap_manifest.py outputs
│   ├── {prefix}_remapped_{assembly}.csv         Full remapped manifest — coordinates + quality columns
│   ├── remapping_Report.txt                     Alignment, pair-filtering, and position-resolution summary
│   └── ambiguous_markers.csv                    Markers with ambiguous mapping (Chr=0)
└── qc/                                          ← qc_filter.py outputs
    ├── matchingSNPs_binary_consistantMapping.{assembly}_map   Main output (final marker map)
    ├── {prefix}_remapped_{assembly}.bim          PLINK BIM format (CHR, SNP, 0, POS, REF, ALT)
    ├── _matchingSNPs_binary_consistantMapping.vcf Final filtered VCF
    ├── _matchingSNPs.vcf                         VCF after design-conflict filter
    ├── _matchingSNPs_binary.vcf                  VCF after polymorphic-site filter
    ├── allele_usage_decision.txt                 Per-SNP orientation decision (as_is / complement)
    ├── QC_Report.txt                             Marker counts at each QC filter stage
    └── remap_assessment/                         MAPQ histograms and known-assembly benchmarks
```

### Remapped CSV Columns

The remapped CSV adds **21 new columns** to every manifest row. All column names that embed the assembly name use the string passed via `-a` (e.g. `-a equCab3` → `Chr_equCab3`).

#### a. Coordinate and position columns

| Column | Type | Meaning |
|---|---|---|
| `Chr_{assembly}` | str | Chromosome (`"0"` = unmapped or ambiguous) |
| `MapInfo_{assembly}` | int | **Final chosen 1-based position** — probe-derived if `CoordDelta < 2`, CIGAR-derived if `CoordDelta ≥ 2` or indel |
| `Strand_{assembly}` | str | `+`, `−`, or `N/A` — TopGenomicSeq alignment strand |
| `Ref_{assembly}` | str | Reference allele in the **TopGenomicSeq alignment orientation**. `qc_filter.py` normalises `Ref` allele to the + strand later for VCF/BIM output. |
| `Alt_{assembly}` | str | Alternate allele in the TopGenomicSeq alignment orientation (same convention as `Ref`) |
| `CoordProbe_{assembly}` | int | Raw probe-derived coordinate before any CIGAR override; `0` for `topseq_only` and unmapped; populated for `probe_only` |
| `Coord_TopSeqCIGAR_{assembly}` | int | CIGAR-walk coordinate from TopGenomicSeq alignment; `0` if SNP in soft clip, `probe_only`, or unmapped |
| `CoordDelta_{assembly}` | int | `\|CoordProbe − Coord_TopSeqCIGAR\|`; `−1` if CIGAR coord unavailable (SNP in soft clip, `topseq_only`, or `probe_only`) |
| `CoordSource_{assembly}` | str | `"probe"` or `"cigar"` — which coordinate is in `MapInfo`; `"N/A"` for unmapped/ambiguous |
| `RefBaseMatch_{assembly}` | str | `"True"` / `"False"` / `"N/A"` — does the genome reference base at `MapInfo` match `Ref` after normalising `Ref` to the + strand? Computed in `remap_manifest.py` as a diagnostic; `qc_filter.py` repeats strand normalisation independently for its design-conflict filter. |
| `ProbeStrand_{assembly}` | str | Alignment strand of the probe: `+`, `−`, or `N/A`. `N/A` for `topseq_only` and unmapped. |
| `StrandAgreementAsExpected_{assembly}` | str | Whether the probe's alignment strand matches the orientation expected from `IlmnStrand`: `"True"`, `"False"`, or `"N/A"`. Always `"True"` or `"N/A"` for `topseq_n_probe` markers (strand disagreement is a hard filter in valid-triple construction). `"N/A"` for rescue paths and unmapped. Used by `qc_filter.py --require-strand-agreement`. |


#### b. Alignment quality columns

| Column | Type | Meaning |
|---|---|---|
| `MAPQ_TopGenomicSeq` | int | MAPQ of winning TopSeq alignment; `NaN` for `probe_only` markers (no TopSeq alignment) |
| `MAPQ_Probe` | int | MAPQ of winning probe alignment; `NaN` for `topseq_only` markers (no probe alignment) |
| `DeltaScore_TopGenomicSeq` | int | AS gap between 1st and 2nd-best TopSeq alignments; `−1` if fewer than 2 alignments |
| `QueryCov_TopGenomicSeq` | float | Fraction of TopSeq query in M/=/X aligned ops (excludes soft/hard clips); `0.0` for unmapped |
| `SoftClipFrac_TopGenomicSeq` | float | Fraction of TopSeq query that is soft-clipped; `0.0` for unmapped |

#### c. Decision columns (see Algorithm Details)

| Column | Values |
|---|---|
| `AlignmentStatus_{assembly}` | `gp1` (both TopSeq alleles + probe), `gp2` (one TopSeq + probe), `gp3` (both TopSeq, no probe), `gp4` (one TopSeq, no probe), `gp5` (probe only, no TopSeq), `unmapped` (nothing aligned). Computed before any filtering. |
| `anchor_{assembly}` | Which source(s) determined the final coordinate: `topseq_n_probe`, `topseq_only`, `probe_only`, `N/A` (see Confidence Tier Summary below) |
| `tie_{assembly}` | How a multi-locus tie was resolved: `unique`, `AS_resolved`, `dAS_resolved`, `NM_resolved`, `CoordDelta_resolved`, `scaffold_resolved`, `ambiguous`, `N/A` |
| `RefAltMethodAgreement_{assembly}` | Agreement between genome ref lookup and NM-based Ref/Alt (see values below) |


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

## Pipeline Decision Tree

The diagram below shows the full per-marker decision flow in `scripts/remap_manifest.py`.

```mermaid
flowchart TD
    classDef tier1    fill:#2d6a4f,color:#fff,stroke:#1b4332,stroke-width:2px
    classDef tier2    fill:#52b788,color:#fff,stroke:#2d6a4f,stroke-width:2px
    classDef tier5    fill:#e63946,color:#fff,stroke:#c1121f,stroke-width:2px
    classDef process  fill:#f8f9fa,color:#000,stroke:#adb5bd,stroke-width:1px

    %% ── Input ──────────────────────────────────────────────────────────────
    INPUT([Illumina Manifest Marker]):::process

    INPUT --> PARSE["extract_candidates\nParse TopGenomicSeq → PREFIX · [A/B] · SUFFIX\nCompute PreLen · PostLen\nDetect assay type I / II\nExpand into TopSeq_A · TopSeq_B sequences"]:::process

    %% ── Alignment ──────────────────────────────────────────────────────────
    PARSE --> ALIGN["minimap2  -ax sr  -N 5\nTopSeq_A + TopSeq_B  (primary + up to 5 secondary)\nAlleleA_ProbeSeq       (primary + up to 5 secondary)\nAll alignments retained — no strand constraint"]:::process

    %% ── Early unmapped exit (nothing aligned at all) ───────────────────────
    ALIGN --> TS_CHECK{"Any alignment\nfound?\n(align_status)"}:::process

    TS_CHECK -->|"unmapped\n(no TopSeq AND no probe\naligned to any locus)"| UNM_TS["unmapped\nChr=0 · Strand=N/A\nno alignment available"]:::tier5

    %% ── Valid-triple construction (gp1–gp5 all proceed here) ───────────────
    TS_CHECK -->|"gp1–gp5\n(≥1 alignment found)"| SBP_BUILD["build_valid_triples + rank_and_resolve\nEnumerate all (TopSeq_allele × TopSeq_align × probe_align) triples\nValidity: TopSeq chr = Probe chr  AND  strand_agreement ∈ {True, N/A}  AND  overlap > 0\nRank: AS_sum → ΔAS_sum → NM_sum → CoordDelta → scaffold → ambiguous\nMAPQ reported as diagnostic only (not used for ranking)"]:::process

    SBP_BUILD --> VALID{"Valid triples\nexist?"}:::process

    %% ── No-valid-triple rescue: TopSeq ─────────────────────────────────────
    VALID -->|"No valid triples\nTopSeq aligned but probe\nabsent · wrong chr · no overlap"| RESCUE["best_topseq_rescue\n_rank_single_aligns waterfall:\nAS → ΔAS → NM → scaffold → ambiguous\ntie: unique · AS_resolved · dAS_resolved\n     NM_resolved · scaffold_resolved · ambiguous"]:::process

    RESCUE --> RESCUE_OUT{"TopSeq rescue\noutcome?"}:::process

    RESCUE_OUT -->|"No TopSeq aligns\nor ambiguous"| PROBE_RESCUE["best_probe_rescue\n_rank_single_aligns waterfall:\nAS → ΔAS → NM → scaffold → ambiguous\ntie: unique · AS_resolved · dAS_resolved\n     NM_resolved · scaffold_resolved · ambiguous"]:::process

    RESCUE_OUT -->|"Winner found\n(placed chr)"| CIGAR_SC["parse_cigar_to_ref_pos\nWalk CIGAR from TopSeq alignment\nQuery index = PreLen (+ strand)\nor PostLen (− strand)"]:::process

    CIGAR_SC --> SC_CHECK{"SNP target\nin soft-clipped\nregion?"}:::process

    SC_CHECK -->|"Yes — target index\nfalls in soft clip\nno ref coord derivable"| UNM_SC["unmapped\nChr=0"]:::tier5

    SC_CHECK -->|"No — CIGAR coord\nsuccessfully derived"| RA_RESCUE["determine_ref_alt_v2\nNM comparison + genome ref lookup\nat rescued position"]:::process

    RA_RESCUE --> RA_RESCUE_CHECK{"Ref/Alt\nresolvable?"}:::process

    RA_RESCUE_CHECK -->|"NM tie\nunresolvable"| UNM_RA["unmapped\nChr=0"]:::tier5

    RA_RESCUE_CHECK -->|"Ref/Alt\nassigned"| TSONLY["topseq_only\nCoord = CIGAR coord\nCoordSource = cigar\nCoordDelta = −1\nMAPQ_Probe = NaN\nanchor = topseq_only"]:::tier2

    %% ── No-valid-triple rescue: Probe ───────────────────────────────────────
    PROBE_RESCUE --> PROBE_OUT{"Probe rescue\noutcome?"}:::process

    PROBE_OUT -->|"No probe\nalignment"| UNM_NO["unmapped\nChr=0"]:::tier5

    PROBE_OUT -->|"Ambiguous"| AMB_PROBE["ambiguous\nChr=0"]:::tier5

    PROBE_OUT -->|"Winner found"| RA_PROBE["determine_ref_alt_v2\nNM comparison + genome ref lookup\nat probe-derived position"]:::process

    RA_PROBE --> PROBEONLY["probe_only\nCoord = probe CIGAR walk\nCoordSource = probe\nMAPQ_TopGenomicSeq = NaN\nanchor = probe_only"]:::tier2

    %% ── Valid-triple resolution ──────────────────────────────────────────────
    VALID -->|"Valid triples\nexist"| RANK["rank_and_resolve\nRanking waterfall:\nAS_sum → ΔAS_sum → NM_sum → CoordDelta → scaffold\ntie values: unique · AS_resolved · dAS_resolved\n            NM_resolved · CoordDelta_resolved · scaffold_resolved\n(supplementary CSV for scaffold_resolved and NM_resolved)"]:::process

    RANK -->|"Winner resolved"| DRA:::process
    RANK -->|"All steps exhausted\ntie = ambiguous\ncompeting rows → ambiguous_markers.csv"| AMB1["ambiguous\nChr=0"]:::tier5

    DRA["determine_ref_alt_v2\nGenome ref lookup (primary) + NM comparison (parallel)\nAgreement reported in RefAltMethodAgreement column"]:::process

    DRA --> NM_WIN{"Ref/Alt\nassignable?"}:::process

    NM_WIN -->|"Yes\n(genome lookup or NM resolved)"| PROBE_COORD["get_probe_coordinate\nProbe CIGAR walk · Probe strand · Assay type\nc_pos = 1-based variant position\n(3′ end for Infinium II · last base for Infinium I)\n(minus-strand deletion correction applied)"]:::process

    NM_WIN -->|"Both methods failed\n(triallelic / no genome match)"| AMB2["ambiguous\nChr=0"]:::tier5

    %% ── CIGAR cross-validation and coordinate selection ─────────────────────
    PROBE_COORD --> CIGAR_XV["parse_cigar_to_ref_pos\nIndependent CIGAR walk on TopSeq alignment\nQuery index = PreLen (+ strand) or PostLen (− strand)\nProduces cigar_coord in parallel with c_pos"]:::process

    CIGAR_XV --> DELTA_CHECK{"CIGAR target\nin soft clip?"}:::process

    DELTA_CHECK -->|"Yes\n(SNP in soft-clipped\nregion of TopSeq)"| USE_PROBE["CoordSource = probe\nCoordDelta = −1\nfinal_pos = c_pos\n(CIGAR unavailable for this marker)"]:::process

    DELTA_CHECK -->|"No\ncigar_coord derived"| DELTA_CALC["CoordDelta = |c_pos − cigar_coord|\nCoordProbe = c_pos  (raw, always stored)\nCoord_TopSeqCIGAR = cigar_coord  (always stored)"]:::process

    DELTA_CALC --> DELTA_THRESH{"CoordDelta\n≥ 2 bp?\nor indel?"}:::process

    DELTA_THRESH -->|"No  (delta = 0 or 1)\nProbe is more reliable\nor both agree"| USE_PROBE2["CoordSource = probe\nfinal_pos = c_pos"]:::process

    DELTA_THRESH -->|"Yes  (delta ≥ 2 bp)\nor indel marker\nCIGAR empirically more accurate"| USE_CIGAR["CoordSource = cigar\nfinal_pos = cigar_coord"]:::process

    USE_PROBE  --> FINAL_COL
    USE_PROBE2 --> FINAL_COL
    USE_CIGAR  --> FINAL_COL

    FINAL_COL["MapInfo_{assembly} = final_pos\nRefBaseMatch = genome base vs strand-normalised Ref\nAll diagnostic columns stored in output CSV"]:::process

    FINAL_COL --> PLACED["mapped / scaffold_resolved\nnm_position_resolved / ref_resolved\nTier 1 — Probe-validated"]:::tier1

    TSONLY --> OUTPUT_TSONLY["topseq_only — Tier 2\nNo probe validation"]:::tier2

    PROBEONLY --> OUTPUT_PROBEONLY["probe_only — Tier 2\nNo TopSeq alignment"]:::tier2

    PLACED  --> OUTPUT_END([Output CSV]):::process
    OUTPUT_TSONLY --> OUTPUT_END
    OUTPUT_PROBEONLY --> OUTPUT_END
```

---

## Algorithm Details

### Infinium Chemistry

- **Infinium I**: two probes (AlleleA_ProbeSeq and AlleleB_ProbeSeq), both ending at the SNP. Variant = **last base** of probe.
- **Infinium II**: `AlleleB_ProbeSeq` is NaN. Variant = base immediately **after** probe 3′ end.
- On the minus strand, the probe's physical 3′ end maps to the alignment **start** position.
- **TopGenomicSeq**: the genomic sequence surrounding the variant position. It has extra IUPAC characters compared to 'SourceSeq' to match the reference genome without indels


### Dual-Alignment Strategy

For each marker, two sequences are aligned to the reference with `minimap2 -ax sr -N 5` (both primary and up to 5 secondary alignments retained):

1. **TopGenomicSeq** — the genomic context `PREFIX[AlleleA/AlleleB]SUFFIX` is split into two candidates (one per allele). The candidate with lower NM at the winning locus is the reference allele.
2. **Probe** (`AlleleA_ProbeSeq`, 50 bp) — aligned independently. Must map to the same chromosome and overlap the TopSeq alignment window; no strand constraint imposed.

A **valid triple** is a `(TopSeq_allele × TopSeq_align × probe_align)` combination that satisfies:
- TopSeq chr == Probe chr
- Strand agreement check (`compute_probe_strand_agreement`) returns `"True"` or `"N/A"`
- Probe–TopSeq overlap > 0 bp

### Tiebreaking Waterfall (main path — `rank_and_resolve`)

When multiple valid triples point to different loci, they are ranked by this cascade until a unique winner is found:

| Step | Criterion | `tie_{assembly}` label |
|---|---|---|
| 1 | All triples at same locus | `unique` |
| 2 | Highest `AS_sum = ts.AS + pb.AS` | `AS_resolved` |
| 3 | Highest `ΔAS_sum` (AS gap vs. competing loci) | `dAS_resolved` |
| 4 | Lowest `NM_sum = ts.NM + pb.NM` | `NM_resolved` |
| 5 | Lowest `CoordDelta` (probe vs. CIGAR coordinate agreement) | `CoordDelta_resolved` |
| 6 | Placed chromosome vs. scaffold | `scaffold_resolved` |
| 7 | All steps exhausted → Chr=0 | `ambiguous` |

> `topseq_only` and `probe_only` markers carry a `tie_{assembly}` value from their single-alignment ranking (AS → ΔAS → NM → scaffold). `CoordDelta_resolved` never appears in rescue paths because there is no probe–CIGAR cross-validation without a valid triple. The per-tie-value breakdown for both rescue paths is in `remapping_Report.txt`.

> **MAPQ is not used for ranking** in the v2 algorithm. It is reported as a diagnostic column only.

### Rescue Paths (no valid triple)

When no valid triple exists, two sequential rescue strategies are attempted:

**TopSeq-only rescue** (`best_topseq_rescue`): uses the same AS → ΔAS → NM → scaffold ranking applied to individual TopSeq alignments (no CoordDelta step — no probe to cross-validate). If a winner is found, the SNP coordinate is derived from a CIGAR walk on the TopSeq alignment. Fails if the SNP target falls in a soft-clipped region or Ref/Alt is unresolvable. Successful markers receive `anchor=topseq_only`, `CoordDelta=−1`, `MAPQ_Probe=NaN`.

**Probe-only rescue** (`best_probe_rescue`): attempted when TopSeq rescue fails or returns ambiguous. Uses the same AS → ΔAS → NM → scaffold ranking on probe alignments. No strand filtering is applied (expected strand cannot be determined without a TopSeq anchor). Successful markers receive `anchor=probe_only`, `CoordDelta=−1`, `MAPQ_TopGenomicSeq=NaN`.



### Coordinate Selection Rule

**`anchor_{assembly}` values:**

| Value | Coordinate evidence |
|---|---|
| `topseq_n_probe` | Probe + TopSeq overlap confirmed; CIGAR cross-check applied (98.7% empirical accuracy) |
| `topseq_only` | No valid triple; TopSeq CIGAR walk used (`CoordSource=cigar`, `MAPQ_Probe=NaN`, `CoordDelta=−1`) (72.9% empirical accuracy) |
| `probe_only` | No TopSeq alignment; probe CIGAR walk used (`CoordSource=probe`, `MAPQ_TopGenomicSeq=NaN`, `CoordDelta=−1`) |
| `N/A` | Chr=0; no reliable genome position assigned |


For `topseq_n_probe` markers, two independent coordinates are computed and compared:

| `CoordDelta` | `MapInfo` source | Rationale |
|---|---|---|
| `0` | Probe | Identical — either works |
| `1` | Probe | Small discrepancy; probe more reliable at delta=1 |
| `≥ 2` | CIGAR | Large discrepancy; CIGAR empirically more accurate (~64–86% vs ~4–17% with probe) |
| `−1` (soft clip) | Probe | CIGAR unavailable; probe used as-is |

Additionally, **indel markers always use the CIGAR coordinate** regardless of `CoordDelta`, because the probe-based coordinate is less precise for insertions/deletions.

`CoordProbe_{assembly}` always stores the raw probe coordinate before any override, enabling retrospective comparison.
`CoordDelta` is a continuous quality signal; `CoordDelta=0` markers are 99.0% coord-accurate.


### Ref/Alt Determination (`determine_ref_alt_v2`)

For each position-resolved marker, two methods run in parallel:

- **Genome lookup** (primary): fetches the reference base at `MapInfo` and compares strand-normalised allele characters against `AlleleA` / `AlleleB`.
- **NM comparison** (parallel): the allele with lower edit distance (NM) at the winning locus is assigned as Ref.

The genome result is used when available. Agreement between the two methods is recorded in `RefAltMethodAgreement_{assembly}`. `NM_unmatch` markers (methods disagree) are worth inspecting for nearby variants.


**`RefAltMethodAgreement_{assembly}` values:**

| Value | Meaning |
|---|---|
| `NM_match` | Genome lookup and NM comparison both succeeded and agree |
| `NM_unmatch` | Both succeeded but disagree — genome result used (flag for QC; inspect for nearby variants) |
| `NM_tied` | Genome succeeded; NM was tied — genome result used |
| `NM_N/A` | Genome succeeded; NM not applicable (`probe_only` marker) |
| `NM_only` | Genome lookup failed; NM result used |
| `ambiguous` | Both methods failed — Chr=0 |

For **indels**: `NM_match` = deletion Ref confirmed by genome fetch; `NM_unmatch` = deletion Ref mismatch (marker removed by design-conflict filter); insertion refs always pass (`NM_match`).

**What is NM?** `NM` is the `NM:i:<n>` edit-distance tag written by minimap2 into each SAM record. It counts mismatches and gap opens between the aligned sequence and the reference — it is **not** derived from CIGAR walking. CIGAR walking is used only for coordinate computation.

---

## QC Filter Cascade

Filters applied sequentially by `qc_filter.py`; `QC_Report.txt` records counts at each stage:

| Stage | Filter condition | Flag |
|---|---|---|
| 1. Unmapped | `Strand_{assembly} == N/A` | always on |
| 2. MAPQ | `MAPQ_TopGenomicSeq < --mapq-topseq` | `--mapq-topseq 30` |
| 2.5. CoordDelta *(optional)* | `CoordDelta > N` OR `anchor_{assembly} == "topseq_only"` | `--coord-delta N` (N ≥ 0); disabled by default |
| 3. Strand agreement *(optional)* | `StrandAgreementAsExpected == False` | `--require-strand-agreement`; disabled by default |
| 4. Design conflict | SNPs: strand-normalised Ref ≠ genome ref base at MapInfo; deletions: pysam fetch of ref ≠ gref; insertions: always pass | always on |
| 5. Exclude indels *(optional)* | Remove all indel markers | `--exclude-indels`; disabled by default |
| 6. Polymorphic sites | Multiple Ref/Alt assignments at the same Chr:Pos | always on |
| 7. Consistency | SAM record count at Chr ≠ 3 (TopSeq_A + TopSeq_B + probe) | requires SAM files |

> The CoordDelta filter explicitly removes `topseq_only` markers (via `anchor_{assembly} == "topseq_only"`) whenever `--coord-delta` is active, because they carry `CoordDelta=−1` and would otherwise numerically pass any threshold ≥ 0.

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

Tests are organised in three files: `tests/test_remap_manifest.py`, `tests/test_qc_filter.py`, and `tests/test_benchmark_compare.py`. Run `pytest --collect-only` to see current totals. Integration tests in `test_benchmark_compare.py` require real data and the `--results-dir` flag.

---

## Benchmarking Remapping Accuracy

For assemblies where the true coordinates are known (e.g., EquCab3-native manifest), run `benchmark_compare.py` after each pipeline execution:

```bash
python scripts/benchmark_compare.py \
    --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
    --assembly  equCab3 \
    --output-dir results_E80selv2_to_equCab3/qc/benchmark/
```

To compare probe-derived, CIGAR-derived, and final coordinates in a three-way accuracy table:

```bash
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
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
