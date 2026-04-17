# Algorithm Overview

How the pipeline decides where each marker goes on the new reference.

---

## The core idea: dual alignment

For every marker, the manifest provides two pieces of DNA:

- **The probe** (`AlleleA_ProbeSeq`, ~50 bp): the actual oligonucleotide that
  hybridises in the assay. High base-level precision but easy to mis-place
  in repetitive regions.
- **The TopGenomicSeq** (~200 bp): the genomic context around the variant
  (`PREFIX[A/B]SUFFIX`). More context → better locus disambiguation, but its
  CIGAR-derived coordinate can be off by a base or two near indels in the
  alignment.

We align **both** to the new reference with `minimap2 -ax sr -N 5` (primary
+ up to 5 secondary alignments). When they agree we have high confidence;
when they disagree the pipeline has well-defined rules for choosing.

---

## Manifest fields the pipeline uses

| Field | Used for |
|---|---|
| `AlleleA_ProbeSeq` | Probe alignment |
| `AlleleB_ProbeSeq` | Distinguishes Infinium I (two probes) from Infinium II (one probe). For Infinium II this is empty. |
| `TopGenomicSeq` | Context alignment; also parsed into `PREFIX`, `[AlleleA / AlleleB]`, `SUFFIX` |
| `IlmnStrand`, `SourceStrand`, `RefStrand` | **Not used by the pipeline.** Manifest design conventions are unreliable as ground truth — see [strand_explained.md](strand_explained.md). `RefStrand` is consumed only by the benchmark for ground-truth strand comparison. |

### Infinium I vs II

- **Infinium I**: two probes (one per allele), both ending **at** the SNP.
  Variant base = the last base of the probe.
- **Infinium II**: one probe ending **just before** the SNP. Variant base =
  the base immediately after the probe's 3′ end.

On the minus strand the probe's 3′ end maps to the alignment **start**
position, not the end.

---

## A "valid triple"

For each marker, we enumerate all combinations of
`(TopSeq_allele × TopSeq_alignment × probe_alignment)`. A combination is
**valid** when:

1. TopSeq and probe align to the same chromosome,
2. The strand-agreement check passes (the probe's alignment strand is
   consistent with the probe-vs-TopSeq orientation evidence — see
   [strand_explained.md](strand_explained.md)),
3. The TopSeq and probe alignment windows overlap (≥ 1 bp).

If at least one valid triple exists, the marker has anchor
`topseq_n_probe`. If not, the pipeline falls back to a rescue path.

---

## The full decision tree

The complete per-marker decision flow is shown as a Mermaid flowchart in
[decision_tree_simple.md](decision_tree_simple.md). Worth reading if you
want to know exactly which branch produced a given output marker.

Briefly, every marker leaves the pipeline labelled by three orthogonal
columns:

- **`anchor_{assembly}`** — which evidence chain placed the marker:
  `topseq_n_probe`, `topseq_only`, `probe_only`, or `N/A` (unmapped).
- **`tie_{assembly}`** — how a multi-locus tie was resolved (or
  `locus_unresolved` if it couldn't be).
- **`RefAltMethodAgreement_{assembly}`** — agreement between the two Ref/Alt
  determination methods (genome lookup vs `NM` comparison).

A marker is fully usable downstream when it has a valid value in all three
of these columns AND `Chr_{assembly}` is not `"0"`.

---

## Coordinate selection (when both probe and TopSeq agree)

For `topseq_n_probe` markers, two coordinates are computed independently:

- **`CoordProbe_{assembly}`** — by walking the **probe's CIGAR** to the variant
  position.
- **`Coord_TopSeqCIGAR_{assembly}`** — by walking the **TopSeq's CIGAR** to the
  variant position.

`CoordDelta_{assembly} = |CoordProbe − Coord_TopSeqCIGAR|`. The chosen final
position is recorded in `MapInfo_{assembly}` and the source in
`CoordSource_{assembly}`:

| `CoordDelta` | Final source (`CoordSource`) | Why |
|---|---|---|
| `0` | `probe_cigar` | Both agree |
| `1` | `probe_cigar` | Small drift; probe slightly more precise at delta = 1 |
| `≥ 2` | `topseq_cigar` | Larger discrepancy; TopSeq CIGAR empirically wins |
| `−1` (CIGAR unavailable, e.g. soft-clipped variant) | `probe_cigar` | TopSeq CIGAR couldn't reach the variant base |

**Indel markers always use `topseq_cigar`** regardless of `CoordDelta` — the
probe coordinate is not precise enough for indels.

For background on CIGAR walking, see [CIGAR_walk.md](CIGAR_walk.md).

---

## Tie-breaking when multiple loci compete

When several valid triples point to different loci, this waterfall picks one:

| Step | Criterion | `tie` label assigned |
|---|---|---|
| 1 | All triples at the same locus | `unique` |
| 2 | Highest `AS_sum = ts.AS + pb.AS` | `AS_resolved` |
| 3 | Highest `ΔAS_sum` (gap to the next-best competitor) | `dAS_resolved` |
| 4 | Lowest `NM_sum = ts.NM + pb.NM` | `NM_resolved` |
| 5 | Lowest `CoordDelta` (probe vs TopSeq agreement) | `CoordDelta_resolved` |
| 6 | Placed chromosome wins over unplaced scaffold | `scaffold_resolved` |
| 7 | All steps exhausted | `locus_unresolved` (Chr=0) |

> MAPQ is reported as a diagnostic column, **not** used for ranking.

The same waterfall (without the CoordDelta step, since there's no probe to
cross-check against) runs in the rescue paths.

---

## Rescue paths

When no valid triple exists:

- **TopSeq-only rescue** (`anchor=topseq_only`) when the TopSeq did align but
  no triple was valid (probe absent, wrong chromosome, no overlap, or
  strand disagreement). Coordinate comes from the TopSeq CIGAR walk.
- **Probe-only rescue** (`anchor=probe_only`) when only the probe aligned (no
  TopSeq alignment at all, "gp5" markers). Coordinate comes from the probe's
  CIGAR walk. No strand-agreement check possible (no TopSeq anchor) so these
  are inherently lower-confidence.
- **TopSeq-rescue ambiguity does not fall through to probe rescue** — once
  TopSeq has failed to resolve the locus, the probe cannot improve on it.

Empirical accuracy on the EquCab2→EquCab3 benchmark:

| Anchor | Markers | Accuracy |
|---|---|---|
| `topseq_n_probe` | majority | ~99.8% |
| `topseq_only` | minority | ~99% |
| `probe_only` | rare | weakest class; segregated in the explanatory benchmark layer |

---

## Ref/Alt determination

For each placed marker, two methods run in parallel:

- **Genome lookup** (primary for SNPs): fetch the reference base at
  `MapInfo`, strand-normalise, compare against `AlleleA` / `AlleleB`.
- **NM comparison** (parallel; primary for indels): the allele with lower
  edit distance (`NM`) at the winning locus is Ref. See
  [NM_comparison.md](NM_comparison.md).

The genome result wins when it's available; otherwise NM is used. The
agreement is recorded in `RefAltMethodAgreement_{assembly}`. The full set
of values, with the variant type each applies to and what the value means:

| Value | Applies to | Meaning |
|---|---|---|
| `NM_match` | SNP | Genome lookup and NM comparison both succeeded and agree. Highest-confidence Ref/Alt outcome. |
| `NM_unmatch` | SNP | Both methods succeeded but disagree on which allele is Ref — genome result is used and recorded; worth inspecting for nearby variants that might have perturbed NM. |
| `NM_tied` | SNP | Genome lookup succeeded; NM comparison produced a tie between alleles — the genome result is used. |
| `NM_N/A` | SNP, Insertion | One method succeeded; the other was unavailable. **SNP:** genome succeeded but no TopSeq alignment was present, so NM couldn't be computed (every `probe_only` marker lands here). **Insertion:** NM assigned Ref/Alt; genome at the variant position was consulted to check for a Ref/Alt swap (see `NM_corrected` below) but it neither confirmed nor contradicted. |
| `NM_only` | SNP | Genome lookup failed (e.g. base wasn't A/C/G/T); the NM result was used as a fallback. |
| `NM_corrected` | Insertion | NM initially assigned Ref = empty string (so the inserted sequence was the Alt), but the genome base at `MapInfo` matched the inserted sequence — meaning the genome has that base, so it is Ref and the missing bases are the variant. Ref/Alt are swapped relative to the NM output. |
| `NM_validated` | Deletion | NM determined which allele is the deleted sequence; the genome was fetched at `MapInfo` (with ±10 bp refinement) and the deletion sequence was confirmed present. High-confidence deletion. |
| `NM_mismatch` | Deletion | NM assigned a deletion sequence but the genome at `MapInfo` does not contain that sequence — the marker is **removed** by Stage 2 of the QC cascade (design conflict filter). |
| `refalt_unresolved` | SNP, Indel | Both methods failed (genome lookup and NM comparison could not determine which allele is Ref). The marker is forced to `Chr=0` and excluded. |
| `N/A` | markers that never reached Ref/Alt determination | The pipeline stopped upstream, so `determine_ref_alt_v2` was never called. Three situations produce this value: (1) the marker was **unmapped** (`anchor == "N/A"`, no valid alignment at all); (2) the marker's **locus was unresolved** (`tie == "locus_unresolved"`, tie-break exhausted in either the main path or a rescue path); (3) the marker's SNP target fell inside a **soft-clipped CIGAR region** in the TopSeq rescue path (no reference coordinate derivable). In all three, `Chr=0`. The 3-D summary table in `remapping_Report.txt` tallies these separately from `refalt_unresolved` (which *is* a Ref/Alt-step failure) — the former appears under `not_attempted(Chr=0)`; the latter under `unresolved(Chr=0)`. |

---

## What can't be remapped

Worked examples of markers the pipeline cannot place — and why those are
limits of the manifest/reference pair, not the algorithm — are in
[cant_remap.md](cant_remap.md). Conversely, markers where the pipeline got it
right and the manifest had it wrong are in [why_we_right.md](why_we_right.md).
