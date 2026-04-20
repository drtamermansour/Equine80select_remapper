# Why We're Right: Markers Where the Manifest Is Wrong

A benchmark-vs-manifest disagreement doesn't always mean the pipeline is wrong.
For a small group of markers out of the non-correct set in the v2 manifest →
EquCab3 run, our pipeline found the true locus and the manifest is the one
that's stale. The worked examples below come from an earlier run that tagged
11 such markers; the individual cases transfer directly even if the exact
count drifts by a few markers when the reference changes.

## How we know

For each non-correct marker we ran a **context check** using the manifest's
`TopGenomicSeq` column (which encodes `PREFIX[A/B]SUFFIX` around the variant):

- `context_forward @ position`: does the 20-bp reference flanking at that position
  equal the manifest's `PREFIX[-20:]` before and `SUFFIX[:20]` after?
- `context_reverse @ position`: same, but against the reverse complement.

A marker is **"we right, manifest wrong"** when context matches at *our*
remapped position but **not** at the *manifest's* position. The 20-bp
flanking is unique enough that a match is strong evidence of the true locus.

## The 11 cases

### Group 1 — Markers whose names literally carry EquCab2 coordinates (4)

These markers were designed against EquCab2 and their `Name` still embeds the
EquCab2 position. The manifest's `Chr` / `MapInfo` columns were never lifted
over to EquCab3, but `TopGenomicSeq` is the canonical DNA context around the
variant. Our remap uses that DNA to find the true EquCab3 position.

| Marker | Manifest (still EquCab2) | Our remap (EquCab3) | Δ |
|---|---|---|---|
| `equCab2:12_22428471_A_C` | chr12:22,428,471 | chr12:25,974,866 | +3.55 Mb |
| `equCab2:12_22428471_A_C_F2BT` | chr12:22,428,471 | chr12:25,897,949 | +3.47 Mb |
| `equCab2:30_12489790_G_A` | chr30:12,489,790 | chr30:12,816,355 | +326 kb |
| `equCab2:30_12489790_G_A_F2BT` | chr30:12,489,790 | chr30:12,935,981 | +446 kb |

Context check for `equCab2:12_22428471_A_C`:
- At chr12:**25,974,866** (ours): the 20-bp flanks match the manifest's
  PREFIX/SUFFIX → forward context verified.
- At chr12:**22,428,471** (manifest): no match in either orientation.

Interpretation: the submitter named the marker after its EquCab2 position
and left the manifest's `Chr`/`MapInfo` columns unchanged. The manifest is
internally inconsistent — `TopGenomicSeq` says one thing, `MapInfo` says
another. Our pipeline trusts the DNA.

### Group 2 — Large-offset silent drift (4)

Non-"equCab2:" markers whose manifest coordinate is also stale, probably from
an older build, with offsets of hundreds of kb to 1+ Mb.

| Marker | Manifest | Our remap | Δ |
|---|---|---|---|
| `Affx-102054887` | chr30:24,214,174 | chr30:25,496,685 | +1.28 Mb |
| `BIEC2_560556` | chr21:29,281,153 | chr21:30,304,193 | +1.02 Mb |
| `BIEC2_821863` | chr30:13,942,656 | chr30:14,388,283 | +446 kb |
| `Affx-101163451` | chr13:8,184,075 | chr13:8,264,674 | +80 kb |

Context forward matches at our position for all four; neither orientation
matches at the manifest position. These were almost certainly correct under
some older reference and were never re-lifted.

### Group 3 — Wrong chromosome in the manifest (2)

Even rarer: the manifest has the marker on the wrong chromosome entirely.

| Marker | Manifest | Our remap |
|---|---|---|
| `CUHSNP00140527` | chr13:19,457,938 | chr10:5,825,043 |
| `BIEC2_1125381` | chrX:59,773,191 | chrX:59,646,167 (reverse context) |

For `CUHSNP00140527`, context forward matches at chr10:5,825,043; no match
anywhere on chr13 near 19,457,938. The manifest's chromosome assignment is
simply wrong.

### Group 4 — Correct locus on an unplaced scaffold (1)

| Marker | Manifest | Our remap |
|---|---|---|
| `BIEC2_709970` | chr26:118,946 | Un_NW_019644702v1:592 (reverse context) |

Context at the scaffold coordinate matches in reverse orientation; nothing at
the chr26 position. The marker's true locus is on an unplaced contig that
represents an alt haplotype not on any placed chromosome.

## What these cases tell us

1. **`TopGenomicSeq` is more reliable than `MapInfo`.** When they disagree,
   trust the sequence. Our whole pipeline is built on that premise and these
   11 cases validate it.
2. **Any manifest carries legacy coordinates.** Names like `equCab2:*` are the
   loudest signal, but even well-behaved names can hide stale positions.
3. **RefStrand can also drift alongside coordinates.** Three of the 11 markers
   disagree with manifest `RefStrand` on top of the coordinate mismatch. The
   benchmark's `manifest_strand_wrong` verdict correctly identifies these
   when forward context matches at our position.

These 11 are the "hero" outputs of the pipeline — cases where running the
remapper produced a *better* answer than what the manifest claimed. They're
counted against us in a naive `correct vs coord_off` metric, but the
context-based explanatory layer pulls them out of that noise.
