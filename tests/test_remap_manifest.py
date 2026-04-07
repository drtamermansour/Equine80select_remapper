import pytest
import sys
import os
import re
import textwrap
from unittest.mock import MagicMock

from remap_manifest import (
    _get_as,
    compute_qcov,
    compute_soft_clip_frac,
    parse_cigar_to_ref_pos,
    is_placed_chromosome,
    _make_competing_rows,
    select_best_pair,
    determine_ref_alt,
    resolve_ref_from_genome,
    parse_topseq_sam,
    get_alignment_end,
    calculate_overlap,
    DecisionCounters,
)


# ── _get_as ───────────────────────────────────────────────────────────────────

def test_get_as_present():
    cols = ["q", "0", "chr1", "100", "60", "50M", "*", "0", "0", "A"*50, "*",
            "NM:i:2", "AS:i:120", "XS:i:80"]
    assert _get_as(cols) == 120

def test_get_as_absent():
    cols = ["q", "0", "chr1", "100", "60", "50M", "*", "0", "0", "A"*50, "*",
            "NM:i:2"]
    assert _get_as(cols) == -1


# ── compute_qcov / compute_soft_clip_frac ─────────────────────────────────────

def test_compute_qcov_full_match():
    """All bases matched → qcov = 1.0."""
    assert compute_qcov("100M") == pytest.approx(1.0)

def test_compute_qcov_with_softclip():
    """10 soft-clipped bases out of 110 total → qcov ≈ 0.909."""
    assert compute_qcov("10S100M") == pytest.approx(100 / 110)

def test_compute_qcov_insertion_excluded_from_numerator():
    """Insertions consume query but are not M/=/X — excluded from qcov numerator."""
    # 90M + 10I + 50M = 150 query bases; only M counts toward aligned → 140/150
    assert compute_qcov("90M10I50M") == pytest.approx(140 / 150)

def test_compute_soft_clip_frac_no_clip():
    assert compute_soft_clip_frac("150M") == pytest.approx(0.0)

def test_compute_soft_clip_frac_both_ends():
    """5S + 140M + 5S = 150 query bases; 10 clipped → 10/150."""
    assert compute_soft_clip_frac("5S140M5S") == pytest.approx(10 / 150)

def test_compute_soft_clip_frac_heavy_clip():
    """50 clipped, 100 aligned → 0.333."""
    assert compute_soft_clip_frac("50S100M") == pytest.approx(50 / 150)


# ── parse_cigar_to_ref_pos ────────────────────────────────────────────────────

def test_cigar_ref_pos_simple_match():
    """All-match CIGAR: query index 10 → POS + 10."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "150M", 10)
    assert pos == 1010
    assert in_sc is False

def test_cigar_ref_pos_leading_softclip_target_after_clip():
    """5S100M: query index 7 (after clip) → POS + (7-5) = POS + 2."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "5S100M", 7)
    assert pos == 1002
    assert in_sc is False

def test_cigar_ref_pos_leading_softclip_target_inside_clip():
    """5S100M: query index 3 (inside soft clip) → in_softclip=True."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "5S100M", 3)
    assert in_sc is True

def test_cigar_ref_pos_deletion_skipped():
    """50M 10D 50M: query index 60 (after deletion) → POS + 50 + 10 + (60-50) = POS + 70."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "50M10D50M", 60)
    assert pos == 1070
    assert in_sc is False

def test_cigar_ref_pos_insertion_target_before():
    """90M 10I 50M: query index 50 (before insertion) → POS + 50."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "90M10I50M", 50)
    assert pos == 1050
    assert in_sc is False

def test_cigar_ref_pos_target_inside_insertion():
    """90M 10I 50M: query index 95 (inside insertion) → POS + 90 (junction), not in_softclip."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "90M10I50M", 95)
    assert pos == 1090
    assert in_sc is False

def test_cigar_ref_pos_minus_strand_postlen_index():
    """Minus-strand: target_idx=PostLen points to first base of RC(allele) in RC query.
    RC query = RC(post)[PostLen] + RC(allele) + RC(pre).
    With PostLen=5, all-match 150M → ref pos POS + 5 = 1005."""
    pos, in_sc = parse_cigar_to_ref_pos(1000, "150M", 5)
    assert pos == 1005
    assert in_sc is False


# ── Alignment dict helpers ────────────────────────────────────────────────────

def _ts(chr, pos, cigar='50M', mapq=60, strand='+', nm=0, as_score=120):
    """Build a minimal TopGenomicSeq alignment dict."""
    return {
        'Chr': chr, 'Pos': pos, 'Cigar': cigar,
        'MAPQ': mapq, 'Strand': strand,
        'End': get_alignment_end(pos, cigar),
        'NM': nm,
        'AS': as_score,
    }


def _pb(chr, pos, cigar='50M', mapq=60, strand='+', as_score=0, nm=0):
    """Build a minimal probe alignment dict."""
    return {
        'Chr': chr, 'Pos': pos, 'Cigar': cigar,
        'MAPQ': mapq, 'Strand': strand,
        'End': get_alignment_end(pos, cigar),
        'AS': as_score,
        'NM': nm,
    }


# ── parse_topseq_sam ──────────────────────────────────────────────────────────

def test_parse_topseq_sam_keeps_secondaries(tmp_path):
    """Secondary alignments (FLAG & 256) must now be retained."""
    sam = textwrap.dedent("""\
        @HD\tVN:1.6
        MarkerA_A\t0\tchr1\t1000\t60\t50M\t*\t0\t0\tACGT\t*\tNM:i:0
        MarkerA_A\t256\tchr5\t2000\t0\t50M\t*\t0\t0\tACGT\t*\tNM:i:5
        MarkerA_B\t0\tchr1\t1000\t55\t50M\t*\t0\t0\tACGT\t*\tNM:i:2
    """)
    sam_file = tmp_path / "test.sam"
    sam_file.write_text(sam)
    result = parse_topseq_sam(str(sam_file))
    assert 'MarkerA' in result
    assert len(result['MarkerA']['A']) == 2       # primary + secondary
    assert len(result['MarkerA']['B']) == 1
    assert result['MarkerA']['A'][0]['Chr'] == 'chr1'
    assert result['MarkerA']['A'][1]['Chr'] == 'chr5'
    assert result['MarkerA']['A'][0]['MAPQ'] == 60
    assert result['MarkerA']['B'][0]['MAPQ'] == 55


def test_parse_topseq_sam_skips_supplementary(tmp_path):
    """Supplementary alignments (FLAG & 2048) must still be skipped."""
    sam = textwrap.dedent("""\
        @HD\tVN:1.6
        MarkerB_A\t0\tchr1\t1000\t60\t50M\t*\t0\t0\tACGT\t*\tNM:i:0
        MarkerB_A\t2048\tchr1\t1200\t60\t30M\t*\t0\t0\tACGT\t*\tNM:i:1
    """)
    sam_file = tmp_path / "test.sam"
    sam_file.write_text(sam)
    result = parse_topseq_sam(str(sam_file))
    assert len(result['MarkerB']['A']) == 1


def test_parse_topseq_sam_parses_as_tag(tmp_path):
    """AS tag must be stored in each alignment dict."""
    sam = textwrap.dedent("""\
        @HD\tVN:1.6
        MarkerX_A\t0\tchr1\t1000\t60\t50M\t*\t0\t0\tACGT\t*\tNM:i:0\tAS:i:145
        MarkerX_B\t0\tchr1\t1000\t55\t50M\t*\t0\t0\tACGT\t*\tNM:i:2
    """)
    sam_file = tmp_path / "test.sam"
    sam_file.write_text(sam)
    result = parse_topseq_sam(str(sam_file))
    assert result['MarkerX']['A'][0]['AS'] == 145
    assert result['MarkerX']['B'][0]['AS'] == -1   # tag absent


# ── is_placed_chromosome ──────────────────────────────────────────────────────

@pytest.mark.parametrize("name", [
    "1", "31", "X", "Y", "MT", "M",
    "chr1", "chr31", "chrX", "chrY", "chrM", "chrMT",
])
def test_is_placed_chromosome_true(name):
    assert is_placed_chromosome(name) is True


@pytest.mark.parametrize("name", [
    "JAAMLG010000001.1",
    "NW_001234567.1",
    "scaffold_42",
    "Un_NW_001234",
    "chrUn_NW_001234",
    "random_contig",
    "",
])
def test_is_placed_chromosome_false(name):
    assert is_placed_chromosome(name) is False


# ── _make_competing_rows ──────────────────────────────────────────────────────

def test_make_competing_rows():
    pairs = [
        ('A', _ts('chr1', 1000, mapq=60, nm=0), _pb('chr1', 1050, mapq=60)),
        ('B', _ts('chr5', 5000, mapq=60, nm=1), _pb('chr5', 5050, mapq=55)),
    ]
    rows = _make_competing_rows(pairs, 'position_tie')
    assert len(rows) == 2
    assert rows[0]['PairRank'] == 1
    assert rows[0]['TopSeqAllele'] == 'A'
    assert rows[0]['AmbiguityReason'] == 'position_tie'
    assert rows[0]['TopSeqChr'] == 'chr1'
    assert rows[0]['TopSeqNM'] == 0
    assert rows[0]['MinMAPQ'] == 60
    assert rows[1]['PairRank'] == 2
    assert rows[1]['MinMAPQ'] == 55


# ── select_best_pair ──────────────────────────────────────────────────────────

def test_select_best_pair_no_valid_pairs_different_chr():
    """Different chromosome between TopSeq and probe → no valid pairs."""
    ts = {'A': [_ts('chr1', 1000, '150M')], 'B': []}
    pb = [_pb('chr2', 1050)]
    assert select_best_pair(ts, pb) is None


def test_select_best_pair_no_valid_pairs_no_overlap():
    """Probe does not overlap TopSeq window → no valid pairs."""
    ts = {'A': [_ts('chr1', 1000, '150M')], 'B': []}
    pb = [_pb('chr1', 2000)]  # far from [1000, 1149]
    assert select_best_pair(ts, pb) is None


def test_select_best_pair_opposite_strand_probe_is_valid():
    """Probe on opposite strand to TopSeq is valid (bottom-strand probe design)."""
    ts = {'A': [_ts('chr1', 1000, '150M', strand='+')], 'B': []}
    pb = [_pb('chr1', 1050, strand='-')]  # overlaps [1000,1149], opposite strand — still valid
    result = select_best_pair(ts, pb)
    assert result is not None
    assert result[0] == 'winner'
    assert result[1] == 'A'


def test_select_best_pair_single_winner():
    """One valid pair → return winner tuple."""
    ts = {'A': [_ts('chr1', 1000, '150M', mapq=60)], 'B': []}
    pb = [_pb('chr1', 1050)]  # overlaps [1000, 1149]
    result = select_best_pair(ts, pb)
    assert result is not None
    assert result[0] == 'winner'
    assert result[1] == 'A'
    assert result[2]['Chr'] == 'chr1'
    assert result[3]['Pos'] == 1050


def test_select_best_pair_redundant_pairs_same_locus():
    """Multiple pairs at the same locus → pick highest combined MAPQ, return winner."""
    ts = {'A': [_ts('chr1', 1000, '150M', mapq=60), _ts('chr1', 1000, '150M', mapq=40)], 'B': []}
    pb = [_pb('chr1', 1050, mapq=60), _pb('chr1', 1060, mapq=40)]
    result = select_best_pair(ts, pb)
    assert result[0] == 'winner'
    assert result[2]['MAPQ'] == 60  # winning topseq
    assert result[3]['MAPQ'] == 60  # winning probe


def test_select_best_pair_scaffold_resolved():
    """One placed chromosome + scaffold at same score → scaffold rule resolves."""
    ts = {
        'A': [
            _ts('chr1', 1000, '150M', mapq=60),
            _ts('JAAMLG010000001.1', 500, '150M', mapq=60),
        ],
        'B': [],
    }
    pb = [
        _pb('chr1', 1050, mapq=60),
        _pb('JAAMLG010000001.1', 550, mapq=60),
    ]
    result = select_best_pair(ts, pb)
    assert result[0] == 'scaffold_resolved'
    assert result[2]['Chr'] == 'chr1'
    assert result[3]['Chr'] == 'chr1'
    competing = result[4]
    assert len(competing) == 2
    assert competing[0]['AmbiguityReason'] == 'scaffold_resolved'


def test_select_best_pair_ambiguous_two_placed_chromosomes():
    """Two placed chromosomes tied → ambiguous."""
    ts = {
        'A': [_ts('chr1', 1000, '150M', mapq=60), _ts('chr5', 5000, '150M', mapq=60)],
        'B': [],
    }
    pb = [_pb('chr1', 1050, mapq=60), _pb('chr5', 5050, mapq=60)]
    result = select_best_pair(ts, pb)
    assert result[0] == 'ambiguous'
    competing = result[1]
    assert len(competing) == 2
    assert competing[0]['AmbiguityReason'] == 'position_tie'


def test_select_best_pair_b_allele_wins():
    """B allele alignment is the only valid pair → winner with allele='B'."""
    ts = {
        'A': [],
        'B': [_ts('chr3', 3000, '150M', mapq=60)],
    }
    pb = [_pb('chr3', 3050, mapq=60)]
    result = select_best_pair(ts, pb)
    assert result[0] == 'winner'
    assert result[1] == 'B'


# ── determine_ref_alt ─────────────────────────────────────────────────────────

def test_determine_ref_alt_a_is_ref():
    """Lower NM on A → A is reference allele."""
    winning_ts = _ts('chr1', 1000, nm=0)
    ts_aligns = {'A': [winning_ts], 'B': [_ts('chr1', 1000, nm=2)]}
    info = {'AlleleA': 'G', 'AlleleB': 'A'}
    result = determine_ref_alt('A', winning_ts, ts_aligns, info)
    assert result == ('G', 'A')


def test_determine_ref_alt_b_is_ref():
    """Lower NM on B → B is reference allele."""
    winning_ts = _ts('chr1', 1000, nm=3)
    ts_aligns = {'A': [winning_ts], 'B': [_ts('chr1', 1000, nm=0)]}
    info = {'AlleleA': 'G', 'AlleleB': 'A'}
    result = determine_ref_alt('A', winning_ts, ts_aligns, info)
    assert result == ('A', 'G')


def test_determine_ref_alt_nm_tie_returns_none():
    """Equal NM between A and B → ambiguous (returns None)."""
    winning_ts = _ts('chr1', 1000, nm=1)
    ts_aligns = {'A': [winning_ts], 'B': [_ts('chr1', 1000, nm=1)]}
    info = {'AlleleA': 'G', 'AlleleB': 'A'}
    assert determine_ref_alt('A', winning_ts, ts_aligns, info) is None


def test_determine_ref_alt_other_allele_absent():
    """Other allele never aligned → winning allele is treated as reference."""
    winning_ts = _ts('chr1', 1000, nm=0)
    ts_aligns = {'A': [winning_ts], 'B': []}
    info = {'AlleleA': 'G', 'AlleleB': 'A'}
    result = determine_ref_alt('A', winning_ts, ts_aligns, info)
    assert result == ('G', 'A')


def test_determine_ref_alt_other_allele_on_different_chr():
    """Other allele aligned only on a different chr → treat as absent."""
    winning_ts = _ts('chr1', 1000, nm=0)
    ts_aligns = {'A': [winning_ts], 'B': [_ts('chr5', 5000, nm=0)]}
    info = {'AlleleA': 'G', 'AlleleB': 'A'}
    result = determine_ref_alt('A', winning_ts, ts_aligns, info)
    assert result == ('G', 'A')


# ── resolve_ref_from_genome ───────────────────────────────────────────────────

def _mock_fasta(return_value):
    """Return a mock pysam.FastaFile whose fetch() returns return_value."""
    fasta = MagicMock()
    fasta.fetch.return_value = return_value
    return fasta


def test_resolve_ref_matches_allele_a():
    """Reference base == allele A → returns (A, B)."""
    fasta = _mock_fasta("A")
    result = resolve_ref_from_genome(fasta, "1", 1000, "A", "G")
    assert result == ("A", "G")
    fasta.fetch.assert_called_once_with("1", 999, 1000)


def test_resolve_ref_matches_allele_b():
    """Reference base == allele B → returns (B, A)."""
    fasta = _mock_fasta("G")
    result = resolve_ref_from_genome(fasta, "5", 500, "A", "G")
    assert result == ("G", "A")


def test_resolve_ref_matches_neither():
    """Reference base == neither allele (true triallelic) → returns None."""
    fasta = _mock_fasta("C")
    result = resolve_ref_from_genome(fasta, "3", 200, "A", "G")
    assert result is None


def test_resolve_ref_fetch_error():
    """pysam raises ValueError on unknown contig → returns None gracefully."""
    fasta = MagicMock()
    fasta.fetch.side_effect = ValueError("unknown contig")
    result = resolve_ref_from_genome(fasta, "chrUn_99", 100, "A", "G")
    assert result is None


# ── DecisionCounters ──────────────────────────────────────────────────────────

def test_decision_counters_ref_resolved_in_summary():
    """format_summary must include ref_resolved line and final_ref_resolved line."""
    c = DecisionCounters()
    c.ref_alt_ref_resolved = 87
    c.final_ref_resolved   = 87
    summary = c.format_summary()
    assert "NM tie resolved by ref lookup" in summary
    assert "87" in summary
    assert "ref-resolved" in summary


# ── select_best_pair NM tiebreaking ──────────────────────────────────────────

def test_select_best_pair_nm_breaks_position_tie():
    """When two placed loci tie on MAPQ, the locus with lower NM wins."""
    ts_aligns = {
        "A": [{"Chr": "1", "Pos": 100, "End": 200, "Strand": "+", "MAPQ": 60, "NM": 0}],
        "B": [{"Chr": "2", "Pos": 300, "End": 400, "Strand": "+", "MAPQ": 60, "NM": 5}],
    }
    pb = [
        {"Chr": "1", "Pos": 150, "End": 200, "Strand": "+", "MAPQ": 60, "Cigar": "50M", "NM": 0},
        {"Chr": "2", "Pos": 350, "End": 400, "Strand": "+", "MAPQ": 60, "Cigar": "50M", "NM": 0},
    ]
    result = select_best_pair(ts_aligns, pb)
    assert result[0] == "nm_position_resolved"
    assert result[2]["Chr"] == "1"   # lower-NM locus wins


def test_select_best_pair_nm_tie_stays_ambiguous():
    """When both placed loci have equal NM, position_tie remains ambiguous."""
    ts_aligns = {
        "A": [{"Chr": "1", "Pos": 100, "End": 200, "Strand": "+", "MAPQ": 60, "NM": 0}],
        "B": [{"Chr": "2", "Pos": 300, "End": 400, "Strand": "+", "MAPQ": 60, "NM": 0}],
    }
    pb = [
        {"Chr": "1", "Pos": 150, "End": 200, "Strand": "+", "MAPQ": 60, "Cigar": "50M", "NM": 0},
        {"Chr": "2", "Pos": 350, "End": 400, "Strand": "+", "MAPQ": 60, "Cigar": "50M", "NM": 0},
    ]
    result = select_best_pair(ts_aligns, pb)
    assert result[0] == "ambiguous"


# ── DecisionCounters nm_position_resolved ────────────────────────────────────

def test_decision_counters_nm_position_resolved_in_summary():
    """format_summary must include nm_position_resolved and nm-position-resolved lines."""
    c = DecisionCounters()
    c.nm_position_resolved = 17
    c.final_nm_position_resolved = 17
    summary = c.format_summary()
    assert "NM position-resolved" in summary
    assert "17" in summary
    assert "nm-position-resolved" in summary


def test_select_best_pair_prefers_higher_as_on_mapq_tie():
    """When min-MAPQ ties, the triple with higher TopSeq AS wins."""
    ts = {
        'A': [
            _ts('chr1', 1000, '150M', mapq=60, nm=0, as_score=150),
            _ts('chr1', 1000, '150M', mapq=60, nm=0, as_score=90),
        ],
        'B': [],
    }
    pb = [_pb('chr1', 1050, mapq=60)]
    result = select_best_pair(ts, pb)
    assert result[0] == 'winner'
    assert result[2]['AS'] == 150
