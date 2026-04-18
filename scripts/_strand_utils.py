"""Strand-normalisation helpers shared across the pipeline.

Single source of truth for the two primitive operations:
  - strand_normalize(allele, strand): rewrite an allele on the + (forward) strand.
    Reverse-complements on '-', passes through on '+' (and any other value).
  - complement(seq): per-base complement (no reversal).

These helpers are imported by remap_manifest.py, qc_filter.py, and
benchmark_compare.py so the same definition is used everywhere.
"""

_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def strand_normalize(allele: str, strand: str) -> str:
    """Return *allele* on the + (forward) strand.

    Reverse-complements when strand == '-'; otherwise returns *allele* unchanged.
    Handles single-base and multi-base alleles uniformly, and preserves case.
    """
    if strand == "-":
        return allele.translate(_RC)[::-1]
    return allele


def complement(seq: str) -> str:
    """Return per-base complement (no reversal)."""
    return seq.translate(_RC)
