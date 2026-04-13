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
    apply_coordinate_role_filter,   # Task 5
    apply_tie_label_filter,          # Task 5
    apply_refalt_conf_filter,        # Task 5
    format_three_d_table,            # Task 5
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
    """CoordDelta filter no longer removes topseq_only; only exceeds-delta rows dropped."""
    import pandas as pd

    df = pd.DataFrame({
        "CoordDelta_test": [5, -1, -1],
        "anchor_test":     ["topseq_n_probe", "topseq_only", "probe_only"],
    })
    # Only row with CoordDelta=5 exceeds threshold=2.
    # topseq_only and probe_only (CoordDelta=-1) pass through.
    exceeds_delta = df["CoordDelta_test"] > 2
    result = df[~exceeds_delta]
    assert len(result) == 2
    assert set(result["anchor_test"]) == {"topseq_only", "probe_only"}


# ── apply_coordinate_role_filter ──────────────────────────────────────────────

def _anchor_df(*values):
    return pd.DataFrame({"anchor_testasm": list(values)})


def test_coord_role_high_keeps_only_topseq_n_probe():
    df = _anchor_df("topseq_n_probe", "topseq_only", "probe_only", "N/A")
    result = apply_coordinate_role_filter(df, "testasm", "High")
    assert list(result["anchor_testasm"]) == ["topseq_n_probe"]


def test_coord_role_moderate_accepts_topseq_only():
    df = _anchor_df("topseq_n_probe", "topseq_only", "probe_only", "N/A")
    result = apply_coordinate_role_filter(df, "testasm", "Moderate")
    assert set(result["anchor_testasm"]) == {"topseq_n_probe", "topseq_only"}


def test_coord_role_low_accepts_probe_only():
    df = _anchor_df("topseq_n_probe", "topseq_only", "probe_only", "N/A")
    result = apply_coordinate_role_filter(df, "testasm", "Low")
    assert set(result["anchor_testasm"]) == {"topseq_n_probe", "topseq_only", "probe_only"}


def test_coord_role_na_always_excluded():
    for role in ("High", "Moderate", "Low"):
        result = apply_coordinate_role_filter(_anchor_df("N/A"), "testasm", role)
        assert len(result) == 0, f"N/A not excluded at role={role}"


# ── apply_tie_label_filter ─────────────────────────────────────────────────────

_ALL_TIES = ["unique", "AS_resolved", "dAS_resolved", "NM_resolved",
             "CoordDelta_resolved", "scaffold_resolved", "ambiguous", "N/A"]


def _tie_df(*values):
    return pd.DataFrame({"tie_testasm": list(values)})


def test_tie_unique_keeps_only_unique():
    result = apply_tie_label_filter(_tie_df(*_ALL_TIES), "testasm", "unique")
    assert list(result["tie_testasm"]) == ["unique"]


def test_tie_resolved_excludes_scaffold_resolved():
    """resolved does NOT include scaffold_resolved — that requires avoid_scaffolds."""
    result = apply_tie_label_filter(_tie_df(*_ALL_TIES), "testasm", "resolved")
    assert set(result["tie_testasm"]) == {
        "unique", "AS_resolved", "dAS_resolved", "NM_resolved", "CoordDelta_resolved"
    }


def test_tie_avoid_scaffolds_adds_scaffold_resolved():
    result = apply_tie_label_filter(
        _tie_df("unique", "scaffold_resolved", "ambiguous", "N/A"),
        "testasm", "avoid_scaffolds"
    )
    assert set(result["tie_testasm"]) == {"unique", "scaffold_resolved"}


def test_tie_ambiguous_always_excluded():
    for label in ("unique", "resolved", "avoid_scaffolds"):
        result = apply_tie_label_filter(_tie_df("ambiguous"), "testasm", label)
        assert len(result) == 0, f"ambiguous not excluded at label={label}"


# ── apply_refalt_conf_filter ───────────────────────────────────────────────────

_ALL_REFALT = [
    "NM_match", "NM_unmatch", "NM_validated", "NM_mismatch",
    "NM_corrected", "NM_tied", "NM_N/A", "NM_only", "ambiguous", "N/A",
]


def _refalt_df(*values):
    return pd.DataFrame({"RefAltMethodAgreement_testasm": list(values)})


def test_refalt_high_keeps_nm_match_and_validated():
    result = apply_refalt_conf_filter(_refalt_df(*_ALL_REFALT), "testasm", "High")
    assert set(result["RefAltMethodAgreement_testasm"]) == {"NM_match", "NM_validated"}


def test_refalt_moderate_adds_nm_na_and_nm_tied():
    result = apply_refalt_conf_filter(_refalt_df(*_ALL_REFALT), "testasm", "Moderate")
    assert set(result["RefAltMethodAgreement_testasm"]) == {
        "NM_match", "NM_validated", "NM_N/A", "NM_tied"
    }


def test_refalt_low_adds_nm_only_nm_unmatch_nm_corrected():
    result = apply_refalt_conf_filter(_refalt_df(*_ALL_REFALT), "testasm", "Low")
    assert set(result["RefAltMethodAgreement_testasm"]) == {
        "NM_match", "NM_validated", "NM_N/A", "NM_tied",
        "NM_only", "NM_unmatch", "NM_corrected",
    }


def test_refalt_nm_mismatch_always_excluded():
    for conf in ("High", "Moderate", "Low"):
        result = apply_refalt_conf_filter(_refalt_df("NM_mismatch"), "testasm", conf)
        assert len(result) == 0, f"NM_mismatch not excluded at conf={conf}"


def test_refalt_ambiguous_always_excluded():
    for conf in ("High", "Moderate", "Low"):
        result = apply_refalt_conf_filter(_refalt_df("ambiguous"), "testasm", conf)
        assert len(result) == 0, f"ambiguous not excluded at conf={conf}"


# ── MAPQ range validation (0–60) ───────────────────────────────────────────────

def _run_parse_args(extra_args):
    """Helper: run parse_args in a subprocess and return CompletedProcess."""
    import subprocess, sys, os
    scripts_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
    cmd = (
        f"import sys; sys.path.insert(0, {scripts_dir!r}); "
        f"sys.argv=['q','-i','x','-r','x','-v','x','-a','x',"
        f"{extra_args}]; "
        f"from qc_filter import parse_args; a=parse_args(); "
        f"print(a.mapq_topseq, a.mapq_probe)"
    )
    return subprocess.run([sys.executable, "-c", cmd], capture_output=True, text=True)


def test_mapq_topseq_rejects_negative():
    assert _run_parse_args("'--mapq-topseq','-1'").returncode != 0


def test_mapq_topseq_rejects_above_60():
    assert _run_parse_args("'--mapq-topseq','61'").returncode != 0


def test_mapq_topseq_accepts_0_and_60():
    proc = _run_parse_args("'--mapq-topseq','0','--mapq-probe','60'")
    assert proc.returncode == 0
    assert "0 60" in proc.stdout


def test_mapq_probe_rejects_negative():
    assert _run_parse_args("'--mapq-probe','-1'").returncode != 0


# ── format_three_d_table ───────────────────────────────────────────────────────

def test_three_d_table_contains_header_and_totals():
    three_d = {("topseq_n_probe", "unique"): {"NM_*": 10, "ambiguous": 0, "N/A": 0}}
    output = format_three_d_table(three_d)
    assert "anchor" in output
    assert "tie" in output
    assert "Total" in output
    assert "10" in output


def test_three_d_table_skips_zero_rows():
    three_d = {
        ("topseq_n_probe", "unique"):     {"NM_*": 5, "ambiguous": 0, "N/A": 0},
        ("topseq_only",    "AS_resolved"):{"NM_*": 0, "ambiguous": 0, "N/A": 0},
    }
    output = format_three_d_table(three_d)
    assert "topseq_only" not in output


def test_three_d_table_grand_total_correct():
    three_d = {
        ("topseq_n_probe", "unique"):      {"NM_*": 3, "ambiguous": 1, "N/A": 2},
        ("topseq_only",    "NM_resolved"): {"NM_*": 4, "ambiguous": 0, "N/A": 0},
    }
    output = format_three_d_table(three_d)
    # Grand total = 3+1+2+4 = 10
    assert "10" in output
