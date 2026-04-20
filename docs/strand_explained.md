# Strand Handling in `remap_manifest.py`

Two completely separate strand concepts are in play, and conflating them was
the source of an earlier ±51 bp bug. The script keeps them strictly separate.

---

## 1. Manifest strand columns — not used by `remap_manifest.py`

| Column | What it encodes | Used by |
|---|---|---|
| `IlmnStrand`   | `TOP` / `BOT` — Illumina design convention | not used by any script |
| `SourceStrand` | `TOP` / `BOT` / `PLUS` / `MINUS` — original design orientation | `benchmark_compare.py` (ground truth only) |
| `RefStrand`    | `+` / `-` | not used |

`remap_manifest.py` does not read any of these three columns. They are
passed through untouched in the output CSV. `SourceStrand` is read later by
`benchmark_compare.py` as strand ground truth (`TOP` / `PLUS` → `+`;
`BOT` / `MINUS` → `-`).

---

## 2. Alignment strand — what `remap_manifest.py` actually uses

Strand information comes entirely from the minimap2 alignments, derived
from the SAM FLAG bit 16:

```python
Strand = "-" if flag & 16 else "+"
```

Two separate alignment strands are tracked per marker:

- `winning_ts["Strand"]` — the TopGenomicSeq alignment strand: `+` or `-`.
- `winning_pb["Strand"]` — the probe alignment strand: `+` or `-`.

These are independent and often opposite — a bottom-strand marker's probe
aligns to the minus strand while its TopGenomicSeq aligns to the plus
strand (or vice versa). There is no strand constraint in pair selection;
a valid pair only requires `chr` match + window overlap.

---

## 3. Where each strand is used and why

### A. `Strand_{assembly}` output column — TopGenomicSeq strand only

```python
new_cols[col_strand].append(winning_ts["Strand"])
```

The strand recorded in the output CSV is always the TopGenomicSeq alignment
strand, never the probe strand. This is the orientation of the genomic
context sequence relative to the reference and is the correct authority
for interpreting `Ref_{assembly}` / `Alt_{assembly}`.

### B. Probe coordinate calculation — probe strand only

```python
c_pos = get_probe_coordinate(
    winning_pb["Pos"], winning_pb["Cigar"], winning_pb["Strand"], assay
)
```

`get_probe_coordinate()` uses exclusively the probe's own strand, never
the TopSeq strand. This was the site of the original ±51 bp bug — the old
code used the TopSeq strand here instead. The probe's strand determines
where its physical 3′ end is:

- `+` strand: 3′ end = `alignment_start + ref_span − 1` (rightmost position)
- `−` strand: 3′ end = `alignment_start` (leftmost position, i.e. POS)

Then, depending on Infinium chemistry:

- Infinium II: variant = base **after** the 3′ end (±1 from probe end)
- Infinium I:  variant = the 3′ end itself

### C. CIGAR coordinate calculation — TopGenomicSeq strand

```python
target_idx = info["PreLen"] if winning_ts["Strand"] == "+" else info["PostLen"]
cigar_coord, cigar_in_sc = parse_cigar_to_ref_pos(
    winning_ts["Pos"], winning_ts["Cigar"], target_idx
)
```

The TopGenomicSeq sequence is `PREFIX[A/B]SUFFIX`. When minimap2 aligns
it to the reference:

- On `+` strand: the sequence is submitted as-is; the SNP bracket starts
  at query index `PreLen`.
- On `−` strand: minimap2 reverse-complements the query internally; the
  bracket now starts at query index `PostLen` (the suffix length becomes
  the leading context in the RC).

The correct query index into the alignment therefore depends on which
strand the TopGenomicSeq aligned to.

### D. Deletion coordinate correction — TopGenomicSeq strand

```python
if len(ref_alt[0]) > len(ref_alt[1]) and winning_ts["Strand"] == "-":
    c_pos -= len(ref_alt[0]) - len(ref_alt[1])
```

For deletions on the minus strand, the probe's 3′ end points to the high
coordinate of the deletion event. Subtracting the deletion length
corrects to the canonical VCF left-anchored position. This only applies
on `−` strand because on `+` strand the probe's 3′ end already points
left (toward the deletion start).

### E. `RefBaseMatch` validation — TopGenomicSeq strand for normalisation

```python
ref_char_fwd = (
    _COMP.get(ref_char, ref_char)
    if winning_ts["Strand"] == "-"
    else ref_char
)
```

`Ref_{assembly}` is stored in alignment strand (the TopGenomicSeq
orientation). To compare it against the forward-strand genome base fetched
via pysam, the script complements it when `Strand == "-"`. This matches
the `strand_normalize()` logic in `qc_filter.py` and correctly predicts
which markers will fail the design-conflict filter.

### F. Expected probe strand — sequence-derived (drives triple validity)

```python
orientation = probe_topseq_orientation(probe_seq, topseq_a, topseq_b)
expected_probe_strand = topseq_strand        if orientation == "same"
                      = flip(topseq_strand)  if orientation == "complement"
agreement = (probe_align_strand == expected_probe_strand)
```

The expected probe-vs-TopSeq strand relationship is derived purely from
sequence comparison — not from `IlmnStrand`. `probe_topseq_orientation`
has two stages:

- **Fast path:** substring presence of `probe_seq` or
  `reverse_complement(probe_seq)` in `topseq_a` / `topseq_b` → `"same"`
  or `"complement"`.
- **Fallback:** 21-mer overlap against `topseq_a` (`_kmer_orientation`)
  when neither substring matches. Ties and degenerate inputs resolve to
  `"same"`.

The fallback guarantees orientation is always resolvable, so `agreement`
is always `"True"` or `"False"` for a winning triple — never `"N/A"`.
`"N/A"` in the `StrandAgreementAsExpected_{assembly}` column only appears
for unmapped / `topseq_only` / `probe_only` markers, which have no probe
or no TopSeq alignment to compare.

Used as a hard filter inside `build_valid_triples`: probes whose
observed alignment strand disagrees with the sequence-derived expectation
are dropped. Reported in the output column
`StrandAgreementAsExpected_{assembly}` for the winning triple
(`topseq_n_probe` anchor).

---

## 4. What `qc_filter.py` does with strand (downstream)

`qc_filter.py` uses `Strand_{assembly}` for strand-normalisation
(`strand_normalize()`) when writing VCF/BIM/map output. The `decision`
column in the final map file (`as_is` / `complement`) is inferred by
matching SNP alleles from the manifest against the strand-normalised
genomic alleles — direct match → `as_is`, complement match → `complement`.
`IlmnStrand`, `SourceStrand`, and `RefStrand` are not consumed by
`qc_filter.py`.

---

## Summary — who owns what

| Strand concept | Source | Used for | Used by |
|---|---|---|---|
| TopGenomicSeq alignment strand (`+`/`-`) | SAM FLAG 16 | output `Strand_{assembly}`, CIGAR target index, deletion correction, `RefBaseMatch`, expected probe strand (combined with sequence orientation) | `remap_manifest.py` + `qc_filter.py` |
| Probe alignment strand (`+`/`-`) | SAM FLAG 16 | `get_probe_coordinate()` and sequence-derived strand-agreement check (hard filter) | `remap_manifest.py` internally |
| `IlmnStrand`   | manifest | nothing (manifest pass-through only)                          | not used |
| `SourceStrand` | manifest | strand ground truth in `benchmark_compare.py` (TOP/PLUS → `+`, BOT/MINUS → `-`) | `benchmark_compare.py` |
| `RefStrand`    | manifest | nothing                                                        | not used |
