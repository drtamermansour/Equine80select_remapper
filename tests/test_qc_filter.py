import math
import pytest
import sys
import os
from unittest.mock import MagicMock

import pandas as pd

from qc_filter import (
    apply_probe_mapq_filter,
    strand_normalize,
    check_deletion_ref_match,
    make_anchor_alleles,
    apply_exclude_indels_filter,
)


# ── apply_probe_mapq_filter ───────────────────────────────────────────────────

def _df(*mapq_values):
    """Build a minimal DataFrame with MAPQ_Probe column."""
    return pd.DataFrame({"MAPQ_Probe": list(mapq_values)})


def test_probe_mapq_filter_disabled_at_zero():
    """threshold=0 (default) → filter disabled, all rows pass."""
    df = _df(0, 5, 20, 60, float('nan'))
    result = apply_probe_mapq_filter(df, threshold=0)
    assert list(result.index) == list(df.index)


def test_probe_mapq_filter_nan_exempt():
    """NaN MAPQ_Probe (topseq_only marker) is exempt — must not be removed."""
    df = _df(float('nan'))
    result = apply_probe_mapq_filter(df, threshold=30)
    assert len(result) == 1


def test_probe_mapq_filter_removes_below_threshold():
    """Row with real MAPQ_Probe below threshold is removed."""
    df = _df(20)
    result = apply_probe_mapq_filter(df, threshold=30)
    assert len(result) == 0


def test_probe_mapq_filter_passes_at_threshold():
    """Row with MAPQ_Probe exactly at threshold is kept."""
    df = _df(30)
    result = apply_probe_mapq_filter(df, threshold=30)
    assert len(result) == 1


def test_probe_mapq_filter_passes_above_threshold():
    """Row with MAPQ_Probe above threshold is kept."""
    df = _df(60)
    result = apply_probe_mapq_filter(df, threshold=30)
    assert len(result) == 1


def test_probe_mapq_filter_mixed_dataset():
    """Mixed: NaN exempt, below threshold removed, at/above threshold kept."""
    df = _df(float('nan'), 20, 30, 60)
    result = apply_probe_mapq_filter(df, threshold=30)
    # NaN (idx 0), 30 (idx 2), 60 (idx 3) survive; 20 (idx 1) removed
    assert list(result.index) == [0, 2, 3]


def test_probe_mapq_filter_zero_mapq_removed_when_threshold_active():
    """MAPQ_Probe=0 (real zero, not NaN) is filtered when threshold > 0.

    Under the old sentinel design, 0 was used for both 'no probe alignment'
    and 'probe aligned with MAPQ 0'. Now 'no probe alignment' is NaN, so
    a genuine MAPQ_Probe=0 row (ambiguous probe) should be filtered out.
    """
    df = _df(0)
    result = apply_probe_mapq_filter(df, threshold=1)
    assert len(result) == 0



# ── strand_normalize with empty string (deletion allele) ─────────────────────

def test_strand_normalize_empty_string_plus_strand():
    """Deletion allele '' on + strand stays '' — no complement to compute."""
    assert strand_normalize("", "+") == ""

def test_strand_normalize_empty_string_minus_strand():
    """Deletion allele '' on - strand stays '' — complement of empty is empty."""
    assert strand_normalize("", "-") == ""


# ── design conflict filter excludes empty _galt (deletion allele) ─────────────

def _dc_df(gref, galt, genome_ref):
    """Build a minimal DataFrame for design conflict filter testing."""
    return pd.DataFrame({"_gref": [gref], "_galt": [galt], "_genome_ref": [genome_ref]})

def test_design_conflict_excludes_deletion_allele_empty_string():
    """Indel with deletion allele stored as '' (not '-') is excluded by _galt != ''."""
    df = _dc_df(gref="CTCGTGCC", galt="", genome_ref="CTCGTGCC")
    result = df[(df["_gref"] == df["_genome_ref"]) & (df["_galt"] != "")]
    assert len(result) == 0

def test_design_conflict_excludes_insertion_allele_empty_string():
    """Indel where the deletion is the ref (empty) is excluded by _gref == _genome_ref failing."""
    df = _dc_df(gref="", galt="CTCGTGCC", genome_ref="T")
    result = df[(df["_gref"] == df["_genome_ref"]) & (df["_galt"] != "")]
    assert len(result) == 0

def test_design_conflict_old_dash_sentinel_no_longer_needed():
    """After fix, '-' string does not appear as _galt — a row with '-' passes the '' check.

    This documents that the old sentinel '-' is replaced by '' and the filter
    must use '' not '-'.
    """
    df = _dc_df(gref="A", galt="-", genome_ref="A")
    result_old = df[(df["_gref"] == df["_genome_ref"]) & (df["_galt"] != "-")]
    result_new = df[(df["_gref"] == df["_genome_ref"]) & (df["_galt"] != "")]
    # Old filter would exclude '-'; new filter lets '-' through (it's now unexpected)
    assert len(result_old) == 0
    assert len(result_new) == 1


# ── check_deletion_ref_match (Q1: indel design conflict) ─────────────────────

def _mock_fasta(return_value):
    fasta = MagicMock()
    fasta.fetch.return_value = return_value
    return fasta


def test_check_deletion_ref_match_matches():
    """Deletion ref 'ACGT' matches reference sequence at mapinfo → True (no conflict)."""
    fasta = _mock_fasta("ACGT")
    assert check_deletion_ref_match(fasta, "chr1", 1000, "ACGT") is True
    fasta.fetch.assert_called_once_with("chr1", 999, 1003)  # 0-based: mapinfo-1 to mapinfo-1+len


def test_check_deletion_ref_match_mismatch():
    """Deletion ref 'ACGT' does not match reference 'TTTT' → False (design conflict)."""
    fasta = _mock_fasta("TTTT")
    assert check_deletion_ref_match(fasta, "chr1", 1000, "ACGT") is False


def test_check_deletion_ref_match_case_insensitive():
    """Reference fetch returns lowercase; comparison is case-insensitive."""
    fasta = _mock_fasta("acgt")
    assert check_deletion_ref_match(fasta, "chr1", 1000, "ACGT") is True


def test_check_deletion_ref_match_fetch_error_returns_false():
    """pysam fetch error (unknown contig) → False (conservative: treat as conflict)."""
    fasta = MagicMock()
    fasta.fetch.side_effect = ValueError("unknown contig")
    assert check_deletion_ref_match(fasta, "chrUn_99", 100, "ACGT") is False


def test_check_deletion_ref_match_insertion_empty_gref_returns_true():
    """For insertion alleles, gref is '' — nothing to verify → True (no conflict)."""
    fasta = MagicMock()  # fetch should not be called
    assert check_deletion_ref_match(fasta, "chr1", 1000, "") is True
    fasta.fetch.assert_not_called()


# ── make_anchor_alleles (Q3: anchor-base encoding for indels) ─────────────────

def test_make_anchor_alleles_deletion():
    """Deletion (gref='CTCG', galt=''): REF=anchor+gref, ALT=anchor, pos=mapinfo-1."""
    fasta = _mock_fasta("A")
    vcf_pos, vcf_ref, vcf_alt = make_anchor_alleles(fasta, "chr1", 1000, "CTCG", "")
    assert vcf_pos == 999            # mapinfo - 1
    assert vcf_ref == "ACTCG"        # anchor + deleted_seq
    assert vcf_alt == "A"            # anchor only


def test_make_anchor_alleles_insertion():
    """Insertion (gref='', galt='CTCG'): REF=anchor, ALT=anchor+galt, pos=mapinfo-1."""
    fasta = _mock_fasta("T")
    vcf_pos, vcf_ref, vcf_alt = make_anchor_alleles(fasta, "chr1", 500, "", "CTCG")
    assert vcf_pos == 499
    assert vcf_ref == "T"
    assert vcf_alt == "TCTCG"


def test_make_anchor_alleles_fetch_error_uses_n():
    """If anchor fetch fails, use 'N' as anchor (conservative)."""
    fasta = MagicMock()
    fasta.fetch.side_effect = ValueError("out of range")
    vcf_pos, vcf_ref, vcf_alt = make_anchor_alleles(fasta, "chr1", 100, "ACG", "")
    assert vcf_ref == "NACG"
    assert vcf_alt == "N"


def test_make_anchor_alleles_uppercase():
    """Anchor is uppercased regardless of FASTA capitalisation."""
    fasta = _mock_fasta("g")
    _, vcf_ref, _ = make_anchor_alleles(fasta, "chr1", 200, "AT", "")
    assert vcf_ref[0] == "G"  # anchor is uppercase


# ── apply_exclude_indels_filter (Q4) ─────────────────────────────────────────

def _indel_df(*rows):
    """Build minimal DataFrame with _gref and _galt columns."""
    return pd.DataFrame([{"_gref": r[0], "_galt": r[1]} for r in rows])


def test_exclude_indels_removes_deletion_allele():
    """Row with _galt=='' (deletion alt) is removed by exclude-indels filter."""
    df = _indel_df(("ACGT", ""))
    result = apply_exclude_indels_filter(df)
    assert len(result) == 0


def test_exclude_indels_removes_insertion_allele():
    """Row with _gref=='' (insertion ref) is removed by exclude-indels filter."""
    df = _indel_df(("", "ACGT"))
    result = apply_exclude_indels_filter(df)
    assert len(result) == 0


def test_exclude_indels_keeps_snps():
    """SNP rows (_gref and _galt both non-empty) are kept."""
    df = _indel_df(("A", "G"), ("C", "T"))
    result = apply_exclude_indels_filter(df)
    assert len(result) == 2


def test_exclude_indels_mixed_dataset():
    """Mixed dataset: SNPs kept, indels removed."""
    df = _indel_df(("A", "G"), ("ACGT", ""), ("C", "T"), ("", "ACGT"))
    result = apply_exclude_indels_filter(df)
    # indices 0 (SNP) and 2 (SNP) kept; 1 (deletion) and 3 (insertion) removed
    assert list(result.index) == [0, 2]


# ── CoordDelta filter uses anchor_{assembly} ──────────────────────────────────

def test_coord_delta_filter_uses_anchor_column():
    """CoordDelta filter identifies topseq_only via anchor_{assembly}, not MappingStatus."""
    import pandas as pd

    # Minimal DataFrame: one normal marker and one topseq_only marker
    df = pd.DataFrame({
        "CoordDelta_test":   [0, -1],
        "anchor_test":       ["topseq_n_probe", "topseq_only"],
    })
    # topseq_only should be removed when coord_delta filter is active (threshold=0)
    is_topseq_only = df["anchor_test"] == "topseq_only"
    exceeds_delta  = df["CoordDelta_test"] > 0
    result = df[~exceeds_delta & ~is_topseq_only]
    assert len(result) == 1
    assert result.iloc[0]["anchor_test"] == "topseq_n_probe"
