# CIGAR Walk: From Query Index to Reference Coordinate

How `parse_cigar_to_ref_pos` in `scripts/remap_manifest.py` turns an
alignment record into a reference coordinate for a specific base in the
query.

## What the walk computes

When minimap2 aligns `seq_A = PREFIX + AlleleA + SUFFIX` to the reference,
the SAM record gives two things:

- **`POS`** (`start_pos`): 1-based reference position of the first
  mapped (non-clipped) query base.
- **`CIGAR`** string: a sequence of operations describing how query bases
  map to reference positions.

The CIGAR walk answers: **"the query base at index `PreLen` (0-based)
landed at what reference coordinate?"**

## CIGAR operations and what they consume

| Op     | Consumes query? | Consumes reference? | Meaning                                 |
|---|---|---|---|
| M/=/X  | yes             | yes                 | Aligned bases (match or mismatch)       |
| I      | yes             | no                  | Query has extra bases not in reference  |
| D / N  | no              | yes                 | Reference has bases the query skipped   |
| S      | yes             | no                  | Soft-clipped (query bases not aligned)  |

The walk maintains two cursors: `curr_q` (query position, 0-based) and
`curr_r` (reference position, 1-based, starts at `POS`). It advances
both according to each operation until `curr_q` reaches
`target_idx = PreLen`.

## Concrete example

```
seq_A = "GATTACA[SNP]GCGTA"     PreLen = 7,  target = seq_A[7] = the SNP base
alignment: POS=1000, CIGAR=5S10M2I8M
```

| Op   | n  | curr_q before | curr_r before | action                                   |
|---|---|---|---|---|
| 5S   | 5  | 0 → 5         | 1000          | `target(7) ≥ 5`, skip                    |
| 10M  | 10 | 5 → 15        | 1000 → 1010   | `target(7) < 15` → **HIT**, return `1000 + (7 − 5) = 1002` |

The SNP base maps to reference position 1002.

## The deletion-allele length-offset fix

For `seq_A = PREFIX + "" + SUFFIX` (`AlleleA` empty), `seq_A[PreLen]` is
`SUFFIX[0]`. For a 43-bp deletion with
`seq_B = PREFIX + CTCGTG...43bp + SUFFIX`:

```
seq_A layout:   [PREFIX][SUFFIX]
                         ↑ PreLen — this is SUFFIX[0]
seq_B layout:   [PREFIX][CTCGTG...43bp][SUFFIX]
                         ↑ PreLen — this is C of CTCGTG
```

The CIGAR walk on `seq_A` returns the reference coordinate of `SUFFIX[0]`.
In the reference genome (which has the full `CTCGTG...43bp`), `SUFFIX[0]`
sits 43 bp past where the deletion sequence starts:

```
Reference:  ...P P P [C T C G T G ... 43bp] [S U F F I X] ...
                      ↑ correct coord         ↑ seq_A lands here (+43)
```

So `cigar_coord = deletion_start + 43`. The fix subtracts
`len(AlleleB) = 43` to recover `deletion_start`.

## The minus-strand subtlety

On the minus strand, minimap2 reports the alignment of `RC(seq_A)`. The
`target_idx` switches to `PostLen` (position of the allele bracket in
the RC sequence), and the returned coordinate is the rightmost reference
base of the allele — not the leftmost. For a single-base SNP this
doesn't matter (left = right), but for multi-base indels on the minus
strand it shifts the coordinate by `allele_len − 1`, which is why the
`CLAUDE.md` pitfall note exists and why
[strand_explained.md § 3.D](strand_explained.md#d-deletion-coordinate-correction--topgenomicseq-strand)
documents the deletion minus-strand correction.

## References in code

- `parse_cigar_to_ref_pos` in `scripts/remap_manifest.py` — the walk itself.
- `get_probe_coordinate` — the probe-side coordinate derivation; uses the
  probe's own strand, unlike the CIGAR walk which uses the TopSeq strand
  (see [strand_explained.md](strand_explained.md)).
