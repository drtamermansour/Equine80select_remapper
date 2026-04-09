import math
import pytest
import sys
import os

import pandas as pd

from qc_filter import apply_probe_mapq_filter, apply_strand_agreement_filter


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


# ── apply_strand_agreement_filter ────────────────────────────────────────────

def _sa_df(*agreement_values):
    """Build a minimal DataFrame with StrandAgreementAsExpected_testasm column."""
    return pd.DataFrame({"StrandAgreementAsExpected_testasm": list(agreement_values)})


def test_strand_agreement_filter_removes_false():
    """StrandAgreementAsExpected=False rows are removed."""
    df = _sa_df("False")
    result = apply_strand_agreement_filter(df, "testasm")
    assert len(result) == 0

def test_strand_agreement_filter_keeps_true():
    """StrandAgreementAsExpected=True rows are kept."""
    df = _sa_df("True")
    result = apply_strand_agreement_filter(df, "testasm")
    assert len(result) == 1

def test_strand_agreement_filter_keeps_na():
    """StrandAgreementAsExpected=N/A rows (topseq_only) are kept — they cannot be evaluated."""
    df = _sa_df("N/A")
    result = apply_strand_agreement_filter(df, "testasm")
    assert len(result) == 1

def test_strand_agreement_filter_mixed():
    """Only False rows removed; True and N/A survive."""
    df = _sa_df("True", "False", "N/A", "True", "False")
    result = apply_strand_agreement_filter(df, "testasm")
    assert list(result.index) == [0, 2, 3]

def test_strand_agreement_filter_wrong_assembly_name_raises():
    """Passing wrong assembly name raises KeyError — guards against silent no-ops."""
    df = _sa_df("True")
    with pytest.raises(KeyError):
        apply_strand_agreement_filter(df, "wrong_asm")
