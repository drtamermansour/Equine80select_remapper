# NM Comparison (Ref/Alt for Indels)

How the pipeline decides which allele is Ref vs Alt for indels, and why the
two internal checks (NM comparison vs deletion validation) can disagree.

See [algorithm_overview.md § Ref/Alt determination](algorithm_overview.md#refalt-determination)
for the higher-level flow and the meaning of each `RefAltMethodAgreement_*`
value.

---

## NM comparison — a relative question

For a deletion marker with `AlleleA = CTCGTG` and `AlleleB = ''`, two
TopSeq sequences are built and aligned:

```
seq_A = PREFIX + CTCGTG + SUFFIX   ← reference allele (has the bases)
seq_B = PREFIX +        + SUFFIX   ← deletion allele (missing the bases)
```

minimap2 aligns both to the EquCab3 reference. The reference genome
contains `CTCGTG` at the wild-type locus (the deletion is a variant, so
the reference is non-deleted). So:

- `seq_A` aligns smoothly → low NM (few edits)
- `seq_B` is missing 6 bases → minimap2 must bridge a gap → higher NM

NM comparison concludes: `AlleleA` is Ref, `AlleleB` is Alt. This is a
**relative judgment**: which allele fits the reference better? It works
even if the coordinate is slightly off, because NM accumulates across
the whole alignment.

---

## Deletion validation — an absolute question

The validation then asks: at the exact derived coordinate `final_pos`,
does `genome[final_pos : final_pos + len(gref)]` equal `gref`?

This is a strict byte-level match at a single position. It fails if:

1. The coordinate is off by even 1 bp (CIGAR walk error), or
2. The EquCab3 sequence at that locus genuinely differs from the
   EquCab2-derived manifest sequence.

---

## Why the two can diverge — a worked example

Take the ACAN D4 marker (`gref =
CTCGTGCCAGATCATCACCACGCAGTCCTCGCCGGCCGTGAAG`, 43 bp):

- **NM comparison:** `seq_A` aligns to EquCab3 chr1 near position
  95257501 with much lower NM than `seq_B` → correctly assigns `AlleleA`
  as Ref.
- **Validation:** fetches 43 bp at chr1:95257501 → gets
  `AAGTTGTCGGGCTGGTTGGGGCGCCAGTTCTCAAATTGCTGTG` → mismatch → marker is
  tagged `NM_mismatch` and removed by QC Stage 2.

The alignment succeeded (minimap2 found a home for the TopSeq near that
locus), NM correctly ranked the alleles, but the exact sequence at the
coordinate doesn't match. This happens because the ACAN gene region in
EquCab3 (Thoroughbred) has diverged from the EquCab2 sequence used to
design the manifest probes. The deletion boundary sequence is different
enough that there's no exact 43-bp match anywhere nearby.

NM tolerated these differences (they span the whole alignment, so
individual base differences don't change which allele wins). The
validation rejects them because it demands an exact match.

---

## Summary

|                                 | NM comparison                                    | Deletion validation                                   |
|---|---|---|
| **Question**                    | Which allele fits better?                        | Does `gref` literally exist at this exact coordinate? |
| **Method**                      | Relative (lower-NM allele wins)                  | Absolute (exact sequence match)                       |
| **Sensitivity to small errors** | Low — distributed across whole alignment         | High — any single mismatch or ±1 bp coord error fails |
| **Can succeed even when…**      | sequence has evolved, coordinate is slightly off | …it cannot                                            |

For small-deletion markers where NM is right and the coordinate is right,
a CIGAR walk that returns position `P + 1` instead of `P` can still make
the validation fail at the wrong spot — the probe coordinate then carries
the correct position.

For divergent structural-variant loci, NM is most likely right about
Ref/Alt, but the EquCab3 sequence has genuinely changed relative to the
EquCab2 manifest design — so the absolute validation cannot pass, and
Stage 2 correctly removes the marker as a design conflict.
