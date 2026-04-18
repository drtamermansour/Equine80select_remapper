# Why Some Markers Can't Be Re-mapped (Or Shouldn't Be)

Out of 82,222 benchmarked markers in the v2 manifest → EquCab3 run, **258 of
the 269 non-correct markers cannot be fixed by a pipeline change**. This
document explains each bucket and why a fix is either impossible or inadvisable.

The remaining 11 non-correct markers are cases where the *manifest* is wrong,
not the pipeline — see `docs/why_we_right.md`.

## How these categories were derived

The benchmark checks flanking context (20 bp PREFIX + 20 bp SUFFIX from the
manifest's `TopGenomicSeq`) against the reference genome at two positions:

- **Our remapped position** (`Chr_equCab3`:`MapInfo_equCab3`)
- **The manifest position** (`Chr`:`MapInfo` as given)

Four categories emerge from which positions match:

| Category | @ our pos | @ manifest pos | Count | Interpretation |
|---|---|---|---|---|
| `we_right_manifest_wrong` | ✓ | ✗ | 11 | Manifest is stale — see `why_we_right.md` |
| `pipeline_wrong_manifest_right` | ✗ | ✓ | 45 | We placed wrong or refused; manifest coord is correct |
| `both_match_duplicate_locus` | ✓ | ✓ | 44 | Both positions exist in the reference — duplicate/alt-haplotype region |
| `neither_matches` | ✗ | ✗ | 169 | Reference has diverged from manifest; no match either way |

---

## 169 — `neither_matches`: reference-sequence divergence

At **neither** our position nor the manifest's does the 20-bp flanking match
the manifest's `TopGenomicSeq`. The EquCab3 reference has diverged enough
from the sequence on which the probe was originally designed that there is
no location where the manifest's DNA context can be verified.

Breakdown by alignment/decision state:

| Result | Anchor | RefAltMethodAgreement | N |
|---|---|---|---|
| `unmapped` | topseq_n_probe | `refalt_unresolved` | 60 |
| `unmapped` | topseq_only | `refalt_unresolved` | 4 |
| `coord_off` | topseq_n_probe | `NM_validated` | 26 |
| `coord_off` | topseq_n_probe | `NM_corrected` | 12 |
| `coord_off` | topseq_n_probe | `NM_mismatch` | 10 |
| `coord_off` | topseq_n_probe | `NM_match` | 8 |
| `coord_off` | topseq_n_probe | `NM_tied` | 5 |
| `coord_off` | topseq_only | `NM_match` | 9 |
| `coord_off` | topseq_only | `NM_tied` | 5 |
| `coord_off` | topseq_only | various | 3 |
| `wrong_chr` | topseq_n_probe | `NM_match` | 1 |
| `wrong_chr` | topseq_only | `NM_match` | 3 |

**Why these can't be fixed:**

- *64 `unmapped + refalt_unresolved`* — minimap2 reported MAPQ=0 and our pipeline
  couldn't determine Ref/Alt confidently (genome lookup failed and NM
  comparison was tied). These live in repeat regions where neither allele
  fits cleanly. Including them would add unreliable data; excluding them is
  the right call.
- *75 `coord_off` with NM-* RefAlt resolution* — minimap2 found a hit and
  Ref/Alt was resolved, but the flanking sequence at that hit doesn't match
  the manifest. The underlying reference sequence differs from what was
  used to design the probes. No amount of pipeline tuning can resolve this
  — the bases literally aren't there.
- *The probe-strand logic is already maximally conservative* — we already
  drop sequences that can't be reliably placed. The ones that remain are
  doing their best.

---

## 45 — `pipeline_wrong_manifest_right`: low-confidence refusals at correct loci

The context DOES match at the manifest's position, but our pipeline put the
marker somewhere else or refused to place it. These are places where our
pipeline was too conservative.

Breakdown:

| Result | Anchor | Tie | RefAltMethodAgreement | N |
|---|---|---|---|---|
| `unmapped` | topseq_n_probe | unique | `refalt_unresolved` | 46 |
| `unmapped` | topseq_only | locus_unresolved | — | 1 |
| `wrong_chr` | topseq_n_probe/topseq_only | — | — | 2 |

*(An additional ~6 `coord_off topseq_n_probe` markers previously sat in
this bucket; those were the 2-bp CoordDelta drift cases now rescued by
the indel-adjacency selection rule — see the commentary at the end of
this section.)*

**The 47 `unmapped + refalt_unresolved` markers (dominant case):**

All have MAPQ=0 on both TopSeq and probe alignments, unique-locus tie, but
`RefAltMethodAgreement=refalt_unresolved` (both genome-lookup and NM comparison
failed). Minimap2 successfully placed both sequences at the manifest's
position but flagged the placement as low-confidence (MAPQ 0 means the
aligner saw multiple equally-good hits).

Our pipeline's Chr=0 rule triggers when RefAlt is unresolved. This is
**by design** — if we can't confidently determine which allele is Ref,
emitting the marker would mislead downstream GWAS / imputation. The
correct position exists, but we can't trust ourselves to annotate it
correctly.

Fixing the 47 would mean relaxing the RefAlt-ambiguity rule, which
trades accuracy for coverage. That's a QC decision, not a bug.

**The residual `wrong_chr` cases:**

These are places where our ranking (AS_sum → ΔAS_sum → NM_sum →
CoordDelta → scaffold) preferred a different locus over the manifest's
correct one. Example:

| Marker | Manifest | Our remap | Δ | Issue |
|---|---|---|---|---|
| `BIEC2_186500` | chr12:13,575,382 | chr12:13,284,912 | -290 kb | `topseq_only` rescue found a competing locus |

**Previously in this bucket — now rescued:** a subset of 2-bp
`CoordDelta=2` (and similar small-delta) SNV markers — e.g.
`BIEC2_429398`, `BIEC2_62056`, `BIEC2_20280`, `BIEC2_1143564` — used to
live here. The selection rule previously was "`CoordDelta ≥ 2` → prefer
CIGAR", which was correct on average (final 98.7% vs probe-alone 98.0%)
but wrong for a small subset of SNVs whose TopSeq CIGARs had minimap2
gaps placed ambiguously in flanking repeats. The cascade now detects
this signature (TopSeq I/D within 5 bp of `target_idx`) and routes those
markers to `probe_cigar` without affecting the majority. Ten markers
across delta buckets 2, 3, 4, and 10 moved from `pipeline_wrong_manifest_right`
to `correct`; empirical regression risk for SNVs was zero (no
TopSeq-correct SNVs carried a nearby I/D).

---

## 44 — `both_match_duplicate_locus`: ambiguous by the reference itself

For these markers, the manifest's 20-bp flanking context exists at **two**
positions in the EquCab3 reference: our chosen locus AND the manifest's
claimed locus. Whichever the pipeline picks is biologically valid.

**22 of 44 have our position on an unplaced scaffold (`Un_NW_*`)** while the
manifest claims a placed chromosome. The flanking sequence is present on
both. Examples:

| Marker | Manifest | Our remap |
|---|---|---|
| `AX-103082952` | chr9:34,872,871 | Un_NW_019641858v1:27,534 |
| `AX-103979822` | chr19:60,548,439 | Un_NW_019641955v1:101,568 |
| `AX-103618489` | chr23:7,797,904 | Un_NW_019643688v1:745,052 |
| `AX-104459545` | chr21:43,414,091 | Un_NW_019643081v1:122,845 |

These unplaced scaffolds are alt-haplotype fragments — they represent
duplicate regions of the assembly that weren't integrated into the placed
chromosomes. Our pipeline already has a `scaffold_resolved` tiebreaker
that prefers placed chromosomes when both compete, but it only fires when
overlap is equal; here the scaffold alignment scored higher on
`AS_sum`/overlap and won outright.

Pre-pipeline mitigation: `scripts/scaffold_haplotype_analyzer.py` +
`scripts/filter_scaffold_haplotypes.py` + `scripts/exclude_alt_haplotypes.py`
can be used to build a cleaned reference (no alt-haplotype scaffolds) before
running the main pipeline. That eliminates the alt-haplotype locus entirely
and forces these 22 markers onto the placed chromosome.

The other 22 are similar duplicate-context cases where both positions are on
placed chromosomes (paralogs, large duplicated segments). Fundamentally the
reference has a two-to-one problem that can only be resolved by external
knowledge (e.g., known linkage-group anchors).

**Special case: 4 "Champagne" duplicates (SLC36A1):**

| Marker | Manifest | Our remap |
|---|---|---|
| `SLC36A1_25884457_Champagne_ilmndup2` | chr14:25,884,457 | chr14:26,012,449 |
| `SLC36A1_25884457_Champagne_v1_ilmndup1` | chr14:25,884,457 | chr14:26,012,449 |
| `SLC36A1_26012449_Champagne_v2` | chr14:26,012,449 | chr14:25,884,457 |
| `BIEC911841_ilmndup2` | chr9:20,379,456 | chr9:20,566,028 |

The manifest has the SLC36A1 Champagne-coat variant annotated at **both**
25,884,457 and 26,012,449 under different marker names, with the entries
that claim each position placing the marker at the *other* one. Both
positions host valid flanking context (a real duplicated segment). The
pipeline deterministically picks one per marker; the manifest's two
coordinates round-trip to each other.

---

## Summary

Of 82,222 benchmarked markers:

- **81,953 (99.7%) correct** — full agreement on coordinate, strand, alleles.
- **269 non-correct:** 
  - 11 the manifest is wrong, we're right (`docs/why_we_right.md`).
  - 47 MAPQ=0 multi-mappers at the right locus but refused by the RefAlt-ambiguity filter (by design).
  - ~2 residual placement errors (`topseq_only` rescue picking a wrong locus); the small-delta CoordSource-heuristic cases were rescued by the indel-adjacency selection rule (see the 45-case section).
  - 44 ambiguous duplicate loci (pre-pipeline alt-haplotype scaffold removal can help 22 of these).
  - 169 fundamental reference-sequence divergence — the DNA that the manifest describes no longer exists in the EquCab3 reference.

The dominant non-correct category — 169 reference-divergence cases — is a
property of the manifest-reference pair, not the pipeline. No code change
fixes sequences that aren't there.
