# Output File Formats

What the pipeline writes, where, and what each column means.

## Directory layout

```
output-dir/
├── temp/                                      intermediate FASTA/SAM (from remap_manifest.py) and _pos/_ref VCFs (from qc_filter.py); removed unless --keep-temp
├── remapping/                                 step 2 — alignment & coordinate resolution
│   ├── {prefix}_remapped_{assembly}.csv         full manifest with quality columns added
│   ├── remapping_Report.txt                     per-decision summary
│   ├── topseq_n_probe_unresolved_markers.csv    competing triples for main-path unresolved markers
│   ├── topseq_only_unresolved_markers.csv       rescue-path unresolved markers (locus or Ref/Alt)
│   ├── probe_only_unresolved_markers.csv        rescue-path unresolved markers (locus or Ref/Alt)
│   ├── scaffold_resolved_markers.csv            competing alignments resolved by placed-chromosome rule
│   └── nm_position_resolved_markers.csv         competing alignments resolved by AS / dAS / NM / CoordDelta
└── qc/                                        step 3 — filter cascade & final outputs
    ├── {prefix}_allele_map_{assembly}.tsv                       ★ main output — manifest↔genome allele crosswalk
    ├── {prefix}_remapped_{assembly}.bim                         PLINK BIM
    ├── {prefix}_remapped_{assembly}.vcf                         final filtered VCF
    ├── {prefix}_remapped_{assembly}_traced.csv                  full manifest + per-marker WhyFiltered column
    ├── QC_Report.txt                                            per-stage filter counts
    └── diagnostics/                                             MAPQ histograms (TopSeq and Probe)
```

---

## Allele map (main output)

`{prefix}_allele_map_{assembly}.tsv` — tab-separated, **with a header row**, 8 columns.
This is the pipeline's headline output: a per-marker crosswalk between the
manifest's original SNP allele encoding and the genomic alleles on the new
reference, with the decision label needed to translate genotype calls from
manifest space to genome space.

| # | Column | Description |
|---|---|---|
| 1 | `chr` | Chromosome (e.g. `1`, `X`, `Un_NW_019641858v1`) |
| 2 | `pos` | 1-based base-pair position |
| 3 | `snp_id` | Marker name from the input manifest |
| 4 | `manifest_alleles` | The two alleles as listed in the manifest's `SNP` column (e.g. `A,G`) |
| 5 | `genomic_alleles` | The same two alleles on the + strand of the reference, in the same order as column 4 |
| 6 | `manifest_ref` | Which of the manifest alleles corresponds to the reference base |
| 7 | `genomic_ref` | The reference base on the + strand |
| 8 | `decision` | How the manifest's allele convention maps to the genome: `as_is` or `complement` (or `indel_as_is` / `indel_complement` for indels) |

Use it alongside a PLINK-format genotype file to convert calls into VCF-standard
Ref/Alt space, or feed the PLINK BIM directly to PLINK2 / Beagle.

---

## Remapped CSV

`{prefix}_remapped_{assembly}.csv` is the input manifest with **21 new columns**
appended. Column names embed the assembly label given via `-a` (e.g. `-a equCab3`
→ `Chr_equCab3`).

### a. Coordinate / position columns

| Column | Type | Meaning |
|---|---|---|
| `Chr_{assembly}` | str | Chromosome (`"0"` = unmapped or unresolved) |
| `MapInfo_{assembly}` | int | **Final 1-based position** chosen by the pipeline |
| `Strand_{assembly}` | str | `+`, `−`, or `N/A` — TopGenomicSeq alignment strand |
| `Ref_{assembly}` | str | Reference allele in the **TopGenomicSeq alignment orientation**. `qc_filter.py` strand-normalises to + strand for VCF/BIM output. |
| `Alt_{assembly}` | str | Alternate allele in the same orientation as `Ref` |
| `CoordProbe_{assembly}` | int | Raw probe-CIGAR coordinate (before any override); `0` if not applicable |
| `Coord_TopSeqCIGAR_{assembly}` | int | TopSeq-CIGAR coordinate; `0` if not applicable |
| `CoordDelta_{assembly}` | float | `\|CoordProbe − Coord_TopSeqCIGAR\|`; `−1` whenever one of the two CIGARs is unavailable — SNP in a soft-clipped TopSeq region, or `topseq_only`, or `probe_only` markers |
| `CoordSource_{assembly}` | str | `"probe_cigar"` (probe alignment's CIGAR) or `"topseq_cigar"` (TopSeq alignment's CIGAR) — which one ended up in `MapInfo`. `"N/A"` for unmapped. |
| `RefBaseMatch_{assembly}` | str | `"True"` / `"False"` / `"N/A"` — does the genome reference base at `MapInfo` match `Ref` after strand normalisation? Diagnostic. |

### b. Alignment quality & diagnostics columns

| Column | Type | Meaning |
|---|---|---|
| `MAPQ_TopGenomicSeq` | int | MAPQ of winning TopSeq alignment; `NaN` for `probe_only` markers |
| `MAPQ_Probe` | int | MAPQ of winning probe alignment; `NaN` for `topseq_only` markers |
| `DeltaScore_TopGenomicSeq` | int | AS gap between best and 2nd-best TopSeq alignments; `−1` if fewer than 2 |
| `QueryCov_TopGenomicSeq` | float | Fraction of TopSeq query in M/=/X aligned ops; `0.0` for unmapped |
| `SoftClipFrac_TopGenomicSeq` | float | Fraction of TopSeq query that is soft-clipped; `0.0` for unmapped |
| `ProbeStrand_{assembly}` | str | Probe alignment strand: `+`, `−`, or `N/A` (`N/A` for `topseq_only` and unmapped) |
| `StrandAgreementAsExpected_{assembly}` | str | `"True"` / `"False"` / `"N/A"` — whether the probe's alignment strand matches what's expected from sequence comparison. Always `"True"` (or `"N/A"` for rescue paths) since the strand check is a hard filter in valid-triple construction. |

### c. Decision columns

| Column | Values |
|---|---|
| `AlignmentStatus_{assembly}` | `gp1` / `gp2` / `gp3` / `gp4` / `gp5` / `unmapped` (raw alignment census; see [algorithm overview](algorithm_overview.md)) |
| `anchor_{assembly}` | `topseq_n_probe` / `topseq_only` / `probe_only` / `N/A` (which evidence chain placed the marker) |
| `tie_{assembly}` | `unique` / `AS_resolved` / `dAS_resolved` / `NM_resolved` / `CoordDelta_resolved` / `scaffold_resolved` / `locus_unresolved` / `N/A` (how multi-locus ties were broken) |
| `RefAltMethodAgreement_{assembly}` | `NM_match` / `NM_validated` / `NM_N/A` / `NM_tied` / `NM_only` / `NM_unmatch` / `NM_corrected` / `NM_mismatch` / `refalt_unresolved` / `N/A` (agreement between genome lookup and NM-based Ref/Alt determination) |

`gp1`–`gp5` are alignment-pattern groups summarising which of the two TopSeq
sequences and/or the probe aligned. Per-group definitions appear in the
"Step 1" section of `remapping_Report.txt`.

For the meaning of each `RefAltMethodAgreement_{assembly}` value, see the
table in [algorithm_overview.md § Ref/Alt determination](algorithm_overview.md#refalt-determination).
For the per-marker decision flow that populates `anchor` / `tie`, see
[algorithm_overview.md](algorithm_overview.md) and
[decision_tree_simple.md](decision_tree_simple.md). Background on the NM
comparison method itself lives in [NM_comparison.md](NM_comparison.md).

---

## Trace CSV

`{prefix}_remapped_{assembly}_traced.csv` is the full input manifest plus one
extra column:

| Column | Meaning |
|---|---|
| `WhyFiltered_{assembly}` | Empty string if the marker passed all filters; otherwise the label of the **first** stage that removed it (e.g. `stage_6_mapq_topseq`, `stage_11_ambiguous_snp`) |

Use this to audit which filter dropped which markers — see
[docs/qc_filters.md](qc_filters.md).

---

## VCF

`{prefix}_remapped_{assembly}.vcf` — the final filtered marker set in VCF v4.3
format. One record per marker surviving all 11 QC stages; chromosome, position,
marker name (ID), REF, and ALT are populated; QUAL/FILTER/INFO are `.`. Indels
use VCF-standard anchor-base encoding (i.e. `pos = mapinfo − 1`,
`REF = anchor + gref`, `ALT = anchor + galt`).

---

## BIM

`{prefix}_remapped_{assembly}.bim` is a standard PLINK BIM file
(chromosome, marker name, cM, position, allele 1, allele 2). Indel rows use
the same anchor-base encoding as the VCF.

---

## QC report

`QC_Report.txt` tabulates how many markers survive after each filter stage,
with the cumulative `(−N)` differences in the right column. The same file
also includes a 3-dimension summary table (anchor × tie × RefAlt outcome) of
the **final** marker set.

---

## Remapping report

`remapping_Report.txt` is the alignment-side equivalent: per-decision
breakdowns for each anchor (`topseq_n_probe`, `topseq_only`, `probe_only`,
`N/A`) and tie-resolution paths, plus diagnostic histograms (MAPQ,
CoordDelta, CoordSource).
