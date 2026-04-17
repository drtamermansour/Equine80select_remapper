import pytest

from _strand_utils import strand_normalize, complement


# ── strand_normalize ─────────────────────────────────────────────────────────

def test_strand_normalize_plus_is_identity_single_base():
    """Plus strand: passed through as-is."""
    assert strand_normalize("A", "+") == "A"


def test_strand_normalize_minus_complements_single_base():
    """Minus strand: single-base RC is single-base complement."""
    assert strand_normalize("A", "-") == "T"
    assert strand_normalize("C", "-") == "G"
    assert strand_normalize("G", "-") == "C"
    assert strand_normalize("T", "-") == "A"


def test_strand_normalize_plus_is_identity_multi_base():
    """Plus strand on a multi-base indel allele: identity."""
    assert strand_normalize("ACGT", "+") == "ACGT"


def test_strand_normalize_minus_reverse_complements_multi_base():
    """Minus strand multi-base: proper reverse-complement (RC), not just complement."""
    assert strand_normalize("ACGT", "-") == "ACGT"      # palindrome
    assert strand_normalize("AAAA", "-") == "TTTT"
    assert strand_normalize("ATCG", "-") == "CGAT"      # complement=TAGC, reverse=CGAT
    assert strand_normalize("CTCGTG", "-") == "CACGAG"  # multi-base deletion case


def test_strand_normalize_empty_allele_insertion_stays_empty():
    """Insertion ref is empty string — must remain empty regardless of strand."""
    assert strand_normalize("", "+") == ""
    assert strand_normalize("", "-") == ""


def test_strand_normalize_lowercase_preserved():
    """Lowercase input: still gets RC-translated, case preserved per str.translate table."""
    assert strand_normalize("acgt", "-") == "acgt"      # lowercase palindrome
    assert strand_normalize("aaaa", "-") == "tttt"


def test_strand_normalize_unknown_strand_is_identity():
    """Defensive: unknown strand value (e.g. 'N/A') passes through without transformation."""
    assert strand_normalize("ACGT", "N/A") == "ACGT"


# ── complement ────────────────────────────────────────────────────────────────

def test_complement_single_base():
    assert complement("A") == "T"
    assert complement("T") == "A"
    assert complement("C") == "G"
    assert complement("G") == "C"


def test_complement_multi_base_no_reverse():
    """complement is complement-only, NOT reverse. Compare to strand_normalize which reverses."""
    assert complement("ATCG") == "TAGC"


def test_complement_empty_string():
    assert complement("") == ""


def test_complement_lowercase():
    assert complement("acgt") == "tgca"
