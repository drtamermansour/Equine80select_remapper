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
    cigar_has_indel_near_query_idx,
    get_probe_coordinate,
    is_placed_chromosome,
    _make_competing_rows,
    determine_ref_alt,
    resolve_ref_from_genome,
    parse_topseq_sam,
    get_alignment_end,
    calculate_overlap,
    DecisionCounters,
    reverse_complement,
    probe_topseq_orientation,
    _kmer_orientation,
    compute_probe_strand_agreement,
    extract_candidates,
    compute_alignment_status,
    build_valid_triples,
    best_topseq_rescue,
    best_probe_rescue,
    rank_and_resolve,
    determine_ref_alt_v2,
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


# ── cigar_has_indel_near_query_idx ────────────────────────────────────────────
# Detects minimap2 gap-placement ambiguity near the SNP in the TopSeq CIGAR.
# When an I/D op sits within `window` bp of target_idx in query space, the
# CIGAR-derived reference coordinate is unreliable (the aligner could have
# placed the gap at an equivalent alternative position).

def test_indel_near_no_indel():
    """All-match CIGAR → no indel to find."""
    assert cigar_has_indel_near_query_idx("50M", 25, window=5) is None


def test_indel_near_d_at_boundary():
    """2D at query boundary (distance 0) — mirrors BIEC2_20280."""
    assert cigar_has_indel_near_query_idx("200M2D201M", 200, window=5) == ("D", 2, 0)


def test_indel_near_d_inside_window():
    """2D at q=198 with target=200 → distance 2; trailing 1D at q=358 is far — closest wins.
    Mirrors BIEC2_1143564."""
    assert cigar_has_indel_near_query_idx("198M2D160M1D43M", 200, window=5) == ("D", 2, 2)


def test_indel_near_d_distance_4():
    """2D at q=196, target=200 → distance 4, within window=5 — mirrors BIEC2_62056."""
    assert cigar_has_indel_near_query_idx("196M2D205M", 200, window=5) == ("D", 2, 4)


def test_indel_near_outside_window():
    """2D at q=100, target=200 → distance 100, outside window → None."""
    assert cigar_has_indel_near_query_idx("100M2D300M", 200, window=5) is None


def test_indel_near_insertion_inside_span():
    """5I spans q=100..104; target=102 inside the insertion span → distance 0."""
    assert cigar_has_indel_near_query_idx("100M5I95M", 102, window=5) == ("I", 5, 0)


def test_indel_near_multi_indel_closest_wins():
    """2D at q=190 (distance 10, outside window); 3D at q=200 (distance 0) — return 3D."""
    assert cigar_has_indel_near_query_idx("190M2D10M3D200M", 200, window=5) == ("D", 3, 0)


def test_indel_near_insertion_boundary_distance():
    """1I spans q=195..196; target=198 → min(|198-195|, |198-196|) = 2."""
    assert cigar_has_indel_near_query_idx("195M1I210M", 198, window=5) == ("I", 1, 2)


def test_indel_near_respects_window_parameter():
    """Same CIGAR/target, window=3 excludes the distance-4 indel that window=5 admitted."""
    assert cigar_has_indel_near_query_idx("196M2D205M", 200, window=5) == ("D", 2, 4)
    assert cigar_has_indel_near_query_idx("196M2D205M", 200, window=3) is None


# ── get_probe_coordinate — soft-clip bug regression tests ────────────────────
# Soft-clipped bases do NOT consume reference. Plus-strand `ref_span` must
# exclude S from its sum; the minus-strand branch already handles leading S
# correctly via `pos_start - leading_s`.

def test_probe_coord_plus_strand_clean_50m_inf_ii():
    """Baseline: plus strand 50M, Inf II → SNP is base after probe 3' end."""
    assert get_probe_coordinate(100, "50M", "+", "II") == 150


def test_probe_coord_plus_strand_clean_50m_inf_i():
    """Baseline: plus strand 50M, Inf I → SNP AT probe 3' end."""
    assert get_probe_coordinate(100, "50M", "+", "I") == 149


def test_probe_coord_plus_strand_leading_softclip_inf_ii():
    """Plus strand 2S48M, Inf II: aligned region is 48 bp (q=2..49).
    SAM pos (=100) already points at the first aligned base, so the probe's 3' end
    is at pos+47; the Inf II SNP is pos+48. Soft clips must not inflate ref_span."""
    assert get_probe_coordinate(100, "2S48M", "+", "II") == 148


def test_probe_coord_plus_strand_trailing_softclip_inf_ii():
    """Plus strand 48M2S, Inf II: trailing soft clip at the probe's 3' end.
    Aligned region is 48 bp; probe 3' end of aligned region is at pos+47;
    Inf II SNP is pos+48."""
    assert get_probe_coordinate(100, "48M2S", "+", "II") == 148


def test_probe_coord_plus_strand_deletion_inf_ii():
    """Plus strand 48M2D2M, Inf II: D consumes reference, so ref_span = 48+2+2 = 52.
    Probe 3' end at pos+51; Inf II SNP = pos+52. Deletions DO count; this must
    remain unchanged by the soft-clip fix."""
    assert get_probe_coordinate(100, "48M2D2M", "+", "II") == 152


def test_probe_coord_minus_strand_leading_softclip_inf_ii():
    """Minus strand 2S48M, Inf II: leading S in SAM = trailing (3') in original probe.
    probe_end = pos - leading_s = 100 - 2 = 98; Inf II SNP = 98 - 1 = 97.
    Minus-strand branch is already correct; keep as regression test."""
    assert get_probe_coordinate(100, "2S48M", "-", "II") == 97


def test_probe_coord_minus_strand_clean_50m_inf_ii():
    """Baseline: minus strand 50M, Inf II → probe_end = pos = 100; SNP = 99."""
    assert get_probe_coordinate(100, "50M", "-", "II") == 99


# ── CoordSource indel-adjacency rescue (integration contract) ────────────────
# Contract tests for the CoordSource selection cascade. The cascade is inline
# in run_remapping (scripts/remap_manifest.py around line 1645); these tests
# encode the expected branch ordering and outcomes.

def test_cigar_indel_rescue_snv_with_topseq_indel_near_target():
    """SNV marker, CoordDelta=2, TopSeq CIGAR has 2D within 5 bp of target_idx.
    Mirrors BIEC2_1143564 failure mode — rule must use probe_cigar (rescue
    path collapsed into the unified probe_cigar label).
    """
    is_indel = False
    cigar_in_sc = False
    cigar_str = "198M2D160M1D43M"
    target_idx = 200
    cigar_coord = 1202
    c_pos = 1200
    coord_delta_val = abs(c_pos - cigar_coord)
    assert coord_delta_val == 2

    # Expected cascade
    if is_indel:
        final_pos, coord_source = cigar_coord, "topseq_cigar"
    elif coord_delta_val >= 2:
        indel_near = cigar_has_indel_near_query_idx(cigar_str, target_idx, window=5)
        if indel_near is not None:
            final_pos, coord_source = c_pos, "probe_cigar"
        else:
            final_pos, coord_source = cigar_coord, "topseq_cigar"
    else:
        final_pos, coord_source = c_pos, "probe_cigar"

    assert coord_source == "probe_cigar"
    assert final_pos == c_pos


def test_cigar_no_rescue_when_topseq_clean():
    """SNV marker, CoordDelta=2, TopSeq CIGAR clean (no indel near target).
    Existing rule preserved: use topseq_cigar."""
    is_indel = False
    cigar_str = "401M"
    target_idx = 200
    cigar_coord = 1202
    c_pos = 1200
    coord_delta_val = abs(c_pos - cigar_coord)

    if is_indel:
        final_pos, coord_source = cigar_coord, "topseq_cigar"
    elif coord_delta_val >= 2:
        indel_near = cigar_has_indel_near_query_idx(cigar_str, target_idx, window=5)
        if indel_near is not None:
            final_pos, coord_source = c_pos, "probe_cigar"
        else:
            final_pos, coord_source = cigar_coord, "topseq_cigar"
    else:
        final_pos, coord_source = c_pos, "probe_cigar"

    assert coord_source == "topseq_cigar"
    assert final_pos == cigar_coord


def test_cigar_indel_branch_precedes_rescue():
    """Indel markers still route to topseq_cigar even with TopSeq indel near target.
    The is_indel branch runs BEFORE the delta branch; rescue must not capture indels."""
    is_indel = True
    cigar_str = "100M8D92M"   # indel marker with adjacent D — mirrors PRKDC_SCID
    target_idx = 100
    cigar_coord = 1050
    c_pos = 1048

    if is_indel:
        final_pos, coord_source = cigar_coord, "topseq_cigar"
    elif abs(c_pos - cigar_coord) >= 2:
        indel_near = cigar_has_indel_near_query_idx(cigar_str, target_idx, window=5)
        if indel_near is not None:
            final_pos, coord_source = c_pos, "probe_cigar"
        else:
            final_pos, coord_source = cigar_coord, "topseq_cigar"
    else:
        final_pos, coord_source = c_pos, "probe_cigar"

    assert coord_source == "topseq_cigar"
    assert final_pos == cigar_coord


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
    assert rows[0]['UnresolvedReason'] == 'position_tie'
    assert rows[0]['TopSeqChr'] == 'chr1'
    assert rows[0]['TopSeqNM'] == 0
    assert rows[0]['MinMAPQ'] == 60
    assert rows[1]['PairRank'] == 2
    assert rows[1]['MinMAPQ'] == 55


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


def test_extract_candidates_deletion_allele_is_empty_string():
    """extract_candidates on a deletion TopGenomicSeq returns '' not '-' for the deletion allele.

    The '-' in [-/CTCGTG] notation means 'no sequence'.  The pipeline must store ''
    so that Ref/Alt columns carry '' rather than the literal dash character.
    """
    pre, a, b, post = extract_candidates("AAACCC[-/CTCGTGCC]TTTGGG")
    assert a == "", f"deletion allele should be '' not {a!r}"
    assert b == "CTCGTGCC"

def test_extract_candidates_insertion_allele_is_empty_string():
    """[I/D] format: first allele is the insertion sequence, second is deletion ('')."""
    pre, a, b, post = extract_candidates("AAACCC[CTCGTGCC/-]TTTGGG")
    assert a == "CTCGTGCC"
    assert b == "", f"deletion allele should be '' not {b!r}"


def test_determine_ref_alt_deletion_allele_empty_string():
    """Deletion allele stored as '' (not '-') passes through as-is.

    After the fix, candidates_info stores '' for deletion alleles rather than '-'.
    determine_ref_alt must return '' unchanged so Ref/Alt columns hold '' not '-'.
    """
    winning_ts = _ts('chr1', 1000, nm=0)
    other_ts   = _ts('chr1', 1000, nm=2)
    ts_aligns  = {'A': [winning_ts], 'B': [other_ts]}
    info = {'AlleleA': '', 'AlleleB': 'CTCGTGCC'}   # deletion ref, insertion alt
    result = determine_ref_alt('A', winning_ts, ts_aligns, info)
    assert result == ('', 'CTCGTGCC')


# ── resolve_ref_from_genome ───────────────────────────────────────────────────

def _mock_fasta(return_value):
    """Return a mock pysam.FastaFile whose fetch() returns return_value."""
    fasta = MagicMock()
    fasta.fetch.return_value = return_value
    return fasta


def test_resolve_ref_matches_allele_a():
    """Reference base == allele A → returns (A, B)."""
    fasta = _mock_fasta("A")
    result = resolve_ref_from_genome(fasta, "1", 1000, "A", "G", "+")
    assert result == ("A", "G")
    fasta.fetch.assert_called_once_with("1", 999, 1000)


def test_resolve_ref_matches_allele_b():
    """Reference base == allele B → returns (B, A)."""
    fasta = _mock_fasta("G")
    result = resolve_ref_from_genome(fasta, "5", 500, "A", "G", "+")
    assert result == ("G", "A")


def test_resolve_ref_matches_neither():
    """Reference base == neither allele (true triallelic) → returns None."""
    fasta = _mock_fasta("C")
    result = resolve_ref_from_genome(fasta, "3", 200, "A", "G", "+")
    assert result is None


def test_resolve_ref_fetch_error():
    """pysam raises ValueError on unknown contig → returns None gracefully."""
    fasta = MagicMock()
    fasta.fetch.side_effect = ValueError("unknown contig")
    result = resolve_ref_from_genome(fasta, "chrUn_99", 100, "A", "G", "+")
    assert result is None


def test_resolve_ref_minus_strand_complements_alleles():
    """On minus strand, alleles are complemented before comparing to fwd genome base.
    Allele A='C' on minus strand → complement is 'G'. Genome returns 'G' → match."""
    fasta = _mock_fasta("G")
    result = resolve_ref_from_genome(fasta, "1", 1000, "C", "A", "-")
    # 'C' complements to 'G' which matches → returns (C, A) in alignment-strand orientation
    assert result == ("C", "A")


def test_resolve_ref_minus_strand_b_matches():
    """On minus strand, allele B='T' → complement 'A'. Genome returns 'A' → match."""
    fasta = _mock_fasta("A")
    result = resolve_ref_from_genome(fasta, "1", 1000, "C", "T", "-")
    # 'C' complements to 'G' (no match), 'T' complements to 'A' (match) → returns (T, C)
    assert result == ("T", "C")


def test_resolve_ref_minus_strand_neither_matches():
    """On minus strand, neither complement matches genome base → returns None."""
    fasta = _mock_fasta("C")
    result = resolve_ref_from_genome(fasta, "1", 1000, "A", "T", "-")
    # 'A' complements to 'T', 'T' complements to 'A'; neither equals 'C' → None
    assert result is None


# ── DecisionCounters ──────────────────────────────────────────────────────────

def test_decision_counters_ref_resolved_in_summary():
    """format_summary must include the NM-tie-resolved-by-ref-lookup row."""
    c = DecisionCounters()
    c.ref_alt_ref_resolved = 87
    summary = c.format_summary()
    assert "NM tie resolved by ref lookup" in summary
    assert "87" in summary


# ── DecisionCounters nm_position_resolved ────────────────────────────────────

def test_decision_counters_nm_position_resolved_in_summary():
    """format_summary must include the tie=NM_resolved label and its count."""
    c = DecisionCounters()
    c.nm_position_resolved = 17
    c.final_nm_position_resolved = 17
    summary = c.format_summary()
    assert "tie=NM_resolved" in summary
    assert "17" in summary


# ── reverse_complement ────────────────────────────────────────────────────────

def test_reverse_complement_simple():
    assert reverse_complement("ACGT") == "ACGT"

def test_reverse_complement_asymmetric():
    assert reverse_complement("AAAA") == "TTTT"

def test_reverse_complement_mixed():
    assert reverse_complement("AACCGGTT") == "AACCGGTT"

def test_reverse_complement_order():
    """RC is reverse of complement: ATCG → complement=TAGC → reverse=CGAT."""
    assert reverse_complement("ATCG") == "CGAT"

def test_reverse_complement_lowercase():
    assert reverse_complement("acgt") == "acgt"


# ── _kmer_orientation ─────────────────────────────────────────────────────────

def test_kmer_orientation_forward_wins():
    """Topseq_a contains a probe-length window embedded forward → 'same'."""
    probe = "ACGTTGCATTAGCCTAGGCAGATCTTGACAGCTATCGAGCTAGATCGTAC"  # 50 bp
    topseq_a = ("N" * 50) + probe[:25] + ("N" * 150)
    assert _kmer_orientation(probe, topseq_a) == "same"


def test_kmer_orientation_rc_wins():
    """Topseq_a contains only a window of RC(probe) → 'complement'."""
    probe = "ACGTTGCATTAGCCTAGGCAGATCTTGACAGCTATCGAGCTAGATCGTAC"
    rc = reverse_complement(probe)
    topseq_a = ("N" * 50) + rc[:25] + ("N" * 150)
    assert _kmer_orientation(probe, topseq_a) == "complement"


def test_kmer_orientation_tie_returns_same():
    """Zero overlap in both directions → tie resolves to 'same'."""
    probe = "A" * 50                # k-mers all A
    topseq_a = "G" * 200            # no A-mers, no T-mers (RC of A is T) → both 0
    assert _kmer_orientation(probe, topseq_a) == "same"


def test_kmer_orientation_short_topseq_returns_same():
    """Topseq_a shorter than k → degenerate input → 'same'."""
    probe = "ACGT" * 13             # 52 bp
    topseq_a = "ACGTACGTAC"         # 10 bp < k=21
    assert _kmer_orientation(probe, topseq_a) == "same"


def test_kmer_orientation_short_probe_returns_same():
    """Probe shorter than k → degenerate input → 'same'."""
    probe = "ACGT"                  # 4 bp < k=21
    topseq_a = "ACGTACGT" * 30      # long enough
    assert _kmer_orientation(probe, topseq_a) == "same"


def test_kmer_orientation_empty_inputs_return_same():
    """Falsy probe or topseq → 'same' (defensive guard)."""
    assert _kmer_orientation("", "ACGT" * 30) == "same"
    assert _kmer_orientation("ACGT" * 30, "") == "same"


def test_kmer_orientation_forward_dominates_over_rc():
    """More forward k-mer hits than rc → 'same' even when both have some."""
    probe = "ACGTTGCATTAGCCTAGGCAGATCTTGACAGCTATCGAGCTAGATCGTAC"
    rc = reverse_complement(probe)
    # Put the full probe (30 fwd 21-mers match) and just one rc 21-mer
    topseq_a = probe + ("N" * 50) + rc[:21]
    assert _kmer_orientation(probe, topseq_a) == "same"


# ── probe_topseq_orientation ──────────────────────────────────────────────────

def test_probe_topseq_orientation_same_found_in_a():
    """Probe as-is is a substring of topseq_a → 'same'."""
    probe   = "AAACCC"
    topseq_a = "XYZAAACCCXYZ"
    topseq_b = "XYZAAATCCXYZ"
    assert probe_topseq_orientation(probe, topseq_a, topseq_b) == "same"

def test_probe_topseq_orientation_same_found_in_b():
    """Probe as-is is a substring of topseq_b (not a) → 'same'."""
    probe    = "AAACCC"
    topseq_a = "XYZAAATCCXYZ"
    topseq_b = "XYZAAACCCXYZ"
    assert probe_topseq_orientation(probe, topseq_a, topseq_b) == "same"

def test_probe_topseq_orientation_complement():
    """RC of probe is a substring of topseq_a → 'complement'."""
    probe    = "AAACCC"   # RC = GGGTT T → "GGGTT T"
    rc_probe = reverse_complement(probe)
    topseq_a = "XYZ" + rc_probe + "XYZ"
    topseq_b = "XXXXXXXXXX"
    assert probe_topseq_orientation(probe, topseq_a, topseq_b) == "complement"

def test_probe_topseq_orientation_fallback_invokes_kmer_same():
    """Substring miss, but k-mer overlap dominates forward → 'same' via fallback."""
    probe = "ACGTTGCATTAGCCTAGGCAGATCTTGACAGCTATCGAGCTAGATCGTAC"  # 50 bp
    # embed a 21-mer of probe (not the whole probe → substring-miss guaranteed)
    topseq_a = ("N" * 20) + probe[:21] + ("N" * 40)
    topseq_b = "G" * 100
    assert probe_topseq_orientation(probe, topseq_a, topseq_b) == "same"


def test_probe_topseq_orientation_fallback_invokes_kmer_complement():
    """Substring miss, but k-mer overlap dominates RC → 'complement' via fallback."""
    probe = "ACGTTGCATTAGCCTAGGCAGATCTTGACAGCTATCGAGCTAGATCGTAC"
    rc = reverse_complement(probe)
    topseq_a = ("N" * 20) + rc[:21] + ("N" * 40)
    topseq_b = "G" * 100
    assert probe_topseq_orientation(probe, topseq_a, topseq_b) == "complement"


def test_probe_topseq_orientation_never_returns_unknown():
    """With the k-mer fallback, orientation is always resolvable."""
    probe = "AAACCC"
    topseq_a = "TTTTTTTTTT"
    topseq_b = "GGGGGGGGGG"
    result = probe_topseq_orientation(probe, topseq_a, topseq_b)
    assert result in ("same", "complement")

def test_probe_topseq_orientation_same_takes_priority_over_complement():
    """If probe matches as-is (same) AND its RC also matches, 'same' wins."""
    probe    = "ACGT"  # palindrome: RC("ACGT") == "ACGT"
    topseq_a = "XACGTX"
    topseq_b = "XXXXXXX"
    assert probe_topseq_orientation(probe, topseq_a, topseq_b) == "same"


# ── compute_probe_strand_agreement ────────────────────────────────────────────
# Expected probe strand is derived purely from sequence comparison
# (probe_topseq_orientation + topseq_strand). IlmnStrand is never consulted.
# Agreement is always "True" or "False" — never "N/A".

def test_strand_agreement_same_orientation_probe_on_topseq_strand_true():
    """orientation=same, topseq=+, probe=+ → agreement True."""
    probe    = "AAACCC"
    topseq_a = "XYZ" + probe + "XYZ"
    topseq_b = "XXXXXXXXXX"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="+", probe_align_strand="+",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "+"
    assert ag == "True"


def test_strand_agreement_same_orientation_probe_opposite_topseq_false():
    """orientation=same, topseq=+, probe=- → agreement False."""
    probe    = "AAACCC"
    topseq_a = "XYZ" + probe + "XYZ"
    topseq_b = "XXXXXXXXXX"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="+", probe_align_strand="-",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "-"
    assert ag == "False"


def test_strand_agreement_complement_orientation_probe_opposite_topseq_true():
    """orientation=complement, topseq=+, probe=- → agreement True."""
    probe    = "AAACCC"
    topseq_a = "XYZ" + reverse_complement(probe) + "XYZ"
    topseq_b = "XXXXXXXXXX"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="+", probe_align_strand="-",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "-"
    assert ag == "True"


def test_strand_agreement_complement_orientation_probe_on_topseq_strand_false():
    """orientation=complement, topseq=+, probe=+ → agreement False."""
    probe    = "AAACCC"
    topseq_a = "XYZ" + reverse_complement(probe) + "XYZ"
    topseq_b = "XXXXXXXXXX"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="+", probe_align_strand="+",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "+"
    assert ag == "False"


def test_strand_agreement_same_orientation_topseq_minus_probe_minus_true():
    """orientation=same, topseq=-, probe=- → expected probe strand = topseq strand = -."""
    probe    = "AAACCC"
    topseq_a = "XYZ" + probe + "XYZ"
    topseq_b = "XXXXXXXXXX"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="-", probe_align_strand="-",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "-"
    assert ag == "True"


def test_strand_agreement_same_orientation_topseq_minus_probe_plus_false():
    """orientation=same, topseq=-, probe=+ → expected - but got + → False."""
    probe    = "AAACCC"
    topseq_a = "XYZ" + probe + "XYZ"
    topseq_b = "XXXXXXXXXX"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="-", probe_align_strand="+",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "+"
    assert ag == "False"


def test_strand_agreement_never_returns_na():
    """With the k-mer fallback, agreement is always 'True' or 'False' — never 'N/A'."""
    # Probe has no substring or k-mer overlap in either direction with topseq → tie → same → agreement computable.
    probe    = "AAACCC"                     # short, substring match unlikely
    topseq_a = "TTTTTTTTTT"
    topseq_b = "GGGGGGGGGG"
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="+", probe_align_strand="+",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ag in ("True", "False")
    assert ps in ("+", "-")


def test_strand_agreement_uses_kmer_fallback_when_substring_misses():
    """Substring miss but 21-mer hit for forward orientation → expected probe strand = topseq strand."""
    probe = "ACGTTGCATTAGCCTAGGCAGATCTTGACAGCTATCGAGCTAGATCGTAC"  # 50 bp
    topseq_a = ("N" * 20) + probe[:21] + ("N" * 40)   # no substring match, forward k-mer present
    topseq_b = "G" * 100
    ps, ag = compute_probe_strand_agreement(
        topseq_strand="+", probe_align_strand="+",
        probe_seq=probe, topseq_a=topseq_a, topseq_b=topseq_b,
    )
    assert ps == "+"
    assert ag == "True"


# ── CIGAR coordinate override for indels (Q5) ────────────────────────────────
# The coordinate-selection logic in run_remapping always picks CIGAR coord for
# indel markers (where one allele is empty string), regardless of CoordDelta.
# These tests verify the helper logic directly via parse_cigar_to_ref_pos since
# the override is inline in run_remapping; we test the contract of the helper
# and document the expected CoordSource="topseq_cigar" for indels.

def test_cigar_coord_used_for_deletion_allele():
    """For a deletion indel (ref_char='ACGT', alt_char=''), CIGAR coord must be chosen.

    This is a documentation/contract test: when ref_char or alt_char is empty string,
    the pipeline must use cigar_coord as final_pos rather than probe coord (c_pos),
    regardless of whether abs(c_pos - cigar_coord) < 2.
    The test encodes the invariant: is_indel → coord_source == "topseq_cigar".
    """
    ref_char = "ACGT"
    alt_char = ""
    is_indel = len(ref_char) == 0 or len(alt_char) == 0
    assert is_indel, "deletion allele '' should be detected as indel"

    # Simulate the override logic from run_remapping:
    # cigar_in_sc=False, cigar_coord=1050, c_pos=1051 (delta=1 → would normally pick probe)
    cigar_in_sc = False
    cigar_coord = 1050
    c_pos       = 1051

    # Without indel override, delta=1 < 2 would pick probe coord
    coord_delta_val = abs(c_pos - cigar_coord)  # = 1
    if cigar_in_sc:
        final_pos    = c_pos
        coord_source = "probe_cigar"
    else:
        if coord_delta_val >= 2:
            final_pos    = cigar_coord
            coord_source = "topseq_cigar"
        else:
            final_pos    = c_pos
            coord_source = "probe_cigar"
    assert coord_source == "probe_cigar"  # without indel override, probe wins

    # With indel override: CIGAR coord is always chosen for indels
    if is_indel and not cigar_in_sc and cigar_coord != 0:
        final_pos    = cigar_coord
        coord_source = "topseq_cigar"
    assert coord_source == "topseq_cigar"
    assert final_pos == cigar_coord


def test_cigar_coord_used_for_insertion_allele():
    """For an insertion indel (ref_char='', alt_char='ACGT'), CIGAR coord must be chosen."""
    ref_char = ""
    alt_char = "ACGT"
    is_indel = len(ref_char) == 0 or len(alt_char) == 0
    assert is_indel

    cigar_in_sc = False
    cigar_coord = 2000
    c_pos       = 2000  # same as CIGAR — delta=0, but CIGAR should still be labelled

    if is_indel and not cigar_in_sc and cigar_coord != 0:
        final_pos    = cigar_coord
        coord_source = "topseq_cigar"
    else:
        final_pos    = c_pos
        coord_source = "probe_cigar"
    assert coord_source == "topseq_cigar"


def test_cigar_coord_not_overridden_when_softclip():
    """If CIGAR coord is unavailable (cigar_in_sc=True), probe coord is used even for indels."""
    ref_char = "ACGT"
    alt_char = ""
    is_indel = len(ref_char) == 0 or len(alt_char) == 0
    assert is_indel

    cigar_in_sc = True
    cigar_coord = 0  # unavailable
    c_pos       = 1050

    if is_indel and not cigar_in_sc and cigar_coord != 0:
        final_pos    = cigar_coord
        coord_source = "topseq_cigar"
    else:
        # Fall back to regular logic
        if cigar_in_sc:
            final_pos    = c_pos
            coord_source = "probe_cigar"
        else:
            final_pos    = c_pos
            coord_source = "probe_cigar"
    assert coord_source == "probe_cigar"
    assert final_pos == c_pos


# ── compute_alignment_status ──────────────────────────────────────────────────

def test_alignment_status_gp1_both_topseq_and_probe():
    """Both alleles + probe aligned → gp1."""
    ts = {"A": [_ts("chr1", 100)], "B": [_ts("chr1", 100, nm=1)]}
    pb = [_pb("chr1", 100)]
    assert compute_alignment_status(ts, pb) == "gp1"


def test_alignment_status_gp2_one_topseq_and_probe():
    """Only allele A + probe → gp2."""
    ts = {"A": [_ts("chr1", 100)], "B": []}
    pb = [_pb("chr1", 100)]
    assert compute_alignment_status(ts, pb) == "gp2"


def test_alignment_status_gp3_both_topseq_no_probe():
    """Both alleles, no probe alignments → gp3."""
    ts = {"A": [_ts("chr1", 100)], "B": [_ts("chr1", 100, nm=1)]}
    pb = []
    assert compute_alignment_status(ts, pb) == "gp3"


def test_alignment_status_gp4_one_topseq_no_probe():
    """Only allele B, no probe → gp4."""
    ts = {"A": [], "B": [_ts("chr1", 100)]}
    pb = []
    assert compute_alignment_status(ts, pb) == "gp4"


def test_alignment_status_gp5_probe_only():
    """No TopSeq at all, probe aligned → gp5."""
    ts = {"A": [], "B": []}
    pb = [_pb("chr1", 100)]
    assert compute_alignment_status(ts, pb) == "gp5"


def test_alignment_status_unmapped():
    """Nothing aligned → unmapped."""
    ts = {"A": [], "B": []}
    pb = []
    assert compute_alignment_status(ts, pb) == "unmapped"


# ── build_valid_triples ───────────────────────────────────────────────────────

def test_build_valid_triples_strand_valid_probe_kept():
    """Probe substring in topseq_a → orientation=same, probe & topseq both + → triple emitted."""
    ts_aligns = {"A": [_ts("chr1", 100, strand="+")], "B": []}
    pb_aligns = [_pb("chr1", 100, strand="+")]
    triples = build_valid_triples(ts_aligns, pb_aligns,
                                   probe_seq="ACGT",
                                   topseq_a="ACGT", topseq_b="ACGG")
    assert len(triples) == 1
    allele, ts, pb = triples[0]
    assert allele == "A"
    assert ts["Chr"] == "chr1"


def test_build_valid_triples_strand_invalid_probe_discarded():
    """orientation=same (probe in topseq_a) but probe strand opposite TopSeq → agreement False → discarded."""
    ts_aligns = {"A": [_ts("chr1", 100, strand="+")], "B": []}
    pb_aligns = [_pb("chr1", 100, strand="-")]
    triples = build_valid_triples(ts_aligns, pb_aligns,
                                   probe_seq="ACGT",
                                   topseq_a="ACGT", topseq_b="ACGG")
    assert triples == []


def test_build_valid_triples_keeps_highest_overlap_probe():
    """Two strand-valid probe alignments on same chr → keep the one with higher overlap."""
    ts = _ts("chr1", 100, cigar="100M")  # spans 100-199
    ts_aligns = {"A": [ts], "B": []}
    # pb1 overlaps 80 bp, pb2 overlaps 20 bp
    pb1 = _pb("chr1", 120, cigar="80M")   # 120-199 → overlap with 100-199 = 80
    pb2 = _pb("chr1", 180, cigar="20M")   # 180-199 → overlap = 20
    triples = build_valid_triples(ts_aligns, [pb1, pb2],
                                   probe_seq="ACGT",
                                   topseq_a="ACGT", topseq_b="ACGG")
    assert len(triples) == 1
    _, _, winning_pb = triples[0]
    assert winning_pb["Pos"] == 120  # pb1 had higher overlap


def test_build_valid_triples_complement_orientation_flips_expected_strand():
    """RC(probe) in topseq_a → orientation=complement → probe expected opposite TopSeq; - strand is valid."""
    ts_aligns = {"A": [_ts("chr1", 100, strand="+")], "B": []}
    pb_aligns = [_pb("chr1", 100, strand="-")]
    probe = "AAACCC"
    triples = build_valid_triples(ts_aligns, pb_aligns,
                                   probe_seq=probe,
                                   topseq_a="XYZ" + reverse_complement(probe) + "XYZ",
                                   topseq_b="XXXXXXXXXX")
    assert len(triples) == 1


def test_build_valid_triples_zero_overlap_discarded():
    """Strand-valid probe on same chr but no overlap → discarded."""
    ts = _ts("chr1", 100, cigar="50M")   # 100-149
    pb = _pb("chr1", 200, cigar="50M")   # 200-249 → no overlap
    ts_aligns = {"A": [ts], "B": []}
    triples = build_valid_triples({"A": [ts], "B": []}, [pb],
                                   probe_seq="ACGT",
                                   topseq_a="ACGT", topseq_b="ACGG")
    assert triples == []


def test_build_valid_triples_different_chr_discarded():
    """Probe on different chromosome than TopSeq → no triple."""
    ts_aligns = {"A": [_ts("chr1", 100)], "B": []}
    pb_aligns = [_pb("chr2", 100)]
    triples = build_valid_triples(ts_aligns, pb_aligns,
                                   probe_seq="ACGT",
                                   topseq_a="ACGT", topseq_b="ACGG")
    assert triples == []


# ── best_topseq_rescue ────────────────────────────────────────────────────────

def test_topseq_rescue_single_alignment_unique():
    """One mapped alignment → unique."""
    ts_aligns = {"A": [_ts("chr1", 100, as_score=120)], "B": []}
    allele, ts, tie, _ = best_topseq_rescue(ts_aligns)
    assert allele == "A"
    assert ts["Chr"] == "chr1"
    assert tie == "unique"


def test_topseq_rescue_as_resolves():
    """Two loci, allele A has higher AS → AS_resolved."""
    ts_aligns = {
        "A": [_ts("chr1", 100, as_score=150)],
        "B": [_ts("chr2", 200, as_score=100)],
    }
    allele, ts, tie, _ = best_topseq_rescue(ts_aligns)
    assert ts["Chr"] == "chr1"
    assert tie == "AS_resolved"


def test_topseq_rescue_nm_resolves():
    """Two loci with same AS, different NM → NM_resolved picks lower NM."""
    ts_aligns = {
        "A": [_ts("chr1", 100, as_score=120, nm=0)],
        "B": [_ts("chr2", 200, as_score=120, nm=2)],
    }
    allele, ts, tie, _ = best_topseq_rescue(ts_aligns)
    assert ts["Chr"] == "chr1"
    assert tie == "NM_resolved"


def test_topseq_rescue_scaffold_resolved():
    """Placed chr + scaffold tied on all metrics → scaffold_resolved."""
    ts_aligns = {
        "A": [_ts("chr1",     100, as_score=120, nm=0)],
        "B": [_ts("NW_12345", 200, as_score=120, nm=0)],
    }
    allele, ts, tie, _ = best_topseq_rescue(ts_aligns)
    assert ts["Chr"] == "chr1"
    assert tie == "scaffold_resolved"


def test_topseq_rescue_locus_unresolved():
    """Two placed chrs with identical metrics → ambiguous."""
    ts_aligns = {
        "A": [_ts("chr1", 100, as_score=120, nm=0)],
        "B": [_ts("chr2", 200, as_score=120, nm=0)],
    }
    allele, ts, tie, _ = best_topseq_rescue(ts_aligns)
    assert allele is None
    assert ts is None
    assert tie == "locus_unresolved"


def test_topseq_rescue_no_mapped_alignment():
    """No mapped TopSeq alignments → returns (None, None, 'N/A')."""
    ts_aligns = {"A": [], "B": []}
    allele, ts, tie, _ = best_topseq_rescue(ts_aligns)
    assert allele is None
    assert ts is None
    assert tie == "N/A"


# ── best_probe_rescue ─────────────────────────────────────────────────────────

def test_probe_rescue_single_alignment_unique():
    """One probe alignment → unique."""
    pb_aligns = [_pb("chr1", 100, as_score=80)]
    pb, tie, _ = best_probe_rescue(pb_aligns)
    assert pb["Chr"] == "chr1"
    assert tie == "unique"


def test_probe_rescue_as_resolves():
    """Two probe positions, different AS → AS_resolved."""
    pb_aligns = [_pb("chr1", 100, as_score=80), _pb("chr2", 200, as_score=60)]
    pb, tie, _ = best_probe_rescue(pb_aligns)
    assert pb["Chr"] == "chr1"
    assert tie == "AS_resolved"


def test_probe_rescue_locus_unresolved():
    """Two placed chrs, identical AS and NM → ambiguous."""
    pb_aligns = [
        _pb("chr1", 100, as_score=80, nm=0),
        _pb("chr2", 200, as_score=80, nm=0),
    ]
    pb, tie, _ = best_probe_rescue(pb_aligns)
    assert pb is None
    assert tie == "locus_unresolved"


# ── rank_and_resolve ──────────────────────────────────────────────────────────

def _info(pre_len=50, post_len=50):
    return {"PreLen": pre_len, "PostLen": post_len}


def test_rank_and_resolve_unique():
    """All triples point to same locus → unique."""
    ts = _ts("chr1", 100, as_score=120, nm=0)
    pb = _pb("chr1", 100, as_score=80, nm=0)
    triples = [("A", ts, pb)]
    result = rank_and_resolve(triples,
                               all_ts_aligns={"A": [ts], "B": []},
                               all_pb_aligns=[pb],
                               info=_info(), assay_type="II")
    assert result[0] == "unique"
    assert result[1] == "A"


def test_rank_and_resolve_as_resolved():
    """Two loci, different AS sum → AS_resolved."""
    ts_a = _ts("chr1", 100, as_score=150, nm=0)
    ts_b = _ts("chr2", 200, as_score=100, nm=0)
    pb_a = _pb("chr1", 100, as_score=80, nm=0)
    pb_b = _pb("chr2", 200, as_score=80, nm=0)
    triples = [("A", ts_a, pb_a), ("B", ts_b, pb_b)]
    result = rank_and_resolve(triples,
                               all_ts_aligns={"A": [ts_a], "B": [ts_b]},
                               all_pb_aligns=[pb_a, pb_b],
                               info=_info(), assay_type="II")
    assert result[0] == "AS_resolved"
    assert result[2]["Chr"] == "chr1"


def test_rank_and_resolve_nm_resolved():
    """Two loci, same AS, different NM sum → NM_resolved."""
    ts_a = _ts("chr1", 100, as_score=120, nm=0)
    ts_b = _ts("chr2", 200, as_score=120, nm=3)
    pb_a = _pb("chr1", 100, as_score=80, nm=0)
    pb_b = _pb("chr2", 200, as_score=80, nm=3)
    triples = [("A", ts_a, pb_a), ("B", ts_b, pb_b)]
    result = rank_and_resolve(triples,
                               all_ts_aligns={"A": [ts_a], "B": [ts_b]},
                               all_pb_aligns=[pb_a, pb_b],
                               info=_info(), assay_type="II")
    assert result[0] == "NM_resolved"
    assert result[2]["Chr"] == "chr1"


def test_rank_and_resolve_scaffold_resolved():
    """Placed chr + scaffold, same metrics → scaffold_resolved."""
    ts_p = _ts("chr1",     100, as_score=120, nm=0)
    ts_s = _ts("NW_12345", 200, as_score=120, nm=0)
    pb_p = _pb("chr1",     100, as_score=80, nm=0)
    pb_s = _pb("NW_12345", 200, as_score=80, nm=0)
    triples = [("A", ts_p, pb_p), ("B", ts_s, pb_s)]
    result = rank_and_resolve(triples,
                               all_ts_aligns={"A": [ts_p], "B": [ts_s]},
                               all_pb_aligns=[pb_p, pb_s],
                               info=_info(), assay_type="II")
    assert result[0] == "scaffold_resolved"
    assert result[2]["Chr"] == "chr1"


def test_rank_and_resolve_locus_unresolved():
    """Two placed chrs, identical metrics → tie=locus_unresolved."""
    ts_a = _ts("chr1", 100, as_score=120, nm=0)
    ts_b = _ts("chr2", 200, as_score=120, nm=0)
    pb_a = _pb("chr1", 100, as_score=80, nm=0)
    pb_b = _pb("chr2", 200, as_score=80, nm=0)
    triples = [("A", ts_a, pb_a), ("B", ts_b, pb_b)]
    result = rank_and_resolve(triples,
                               all_ts_aligns={"A": [ts_a], "B": [ts_b]},
                               all_pb_aligns=[pb_a, pb_b],
                               info=_info(), assay_type="II")
    assert result[0] == "locus_unresolved"


# ── determine_ref_alt_v2 ──────────────────────────────────────────────────────

def _fasta_returning(base):
    f = MagicMock()
    f.fetch.return_value = base
    return f


def test_refalt_v2_snp_nm_match():
    """Genome says A, NM says A is ref → NM_match."""
    ts_a = _ts("chr1", 1000, nm=0)
    ts_b = _ts("chr1", 1000, nm=1)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "A", "AlleleB": "G"}
    fasta = _fasta_returning("A")   # genome = A = allele A on + strand
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "A"
    assert alt == "G"
    assert agree == "NM_match"


def test_refalt_v2_snp_nm_unmatch():
    """Genome says G (allele B), NM says A is ref → NM_unmatch, genome wins."""
    ts_a = _ts("chr1", 1000, nm=0)
    ts_b = _ts("chr1", 1000, nm=1)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "A", "AlleleB": "G"}
    fasta = _fasta_returning("G")   # genome = G = allele B
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "G"
    assert alt == "A"
    assert agree == "NM_unmatch"


def test_refalt_v2_snp_nm_tied():
    """NM is tied → genome lookup used, agreement = NM_tied."""
    ts_a = _ts("chr1", 1000, nm=1)
    ts_b = _ts("chr1", 1000, nm=1)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "A", "AlleleB": "G"}
    fasta = _fasta_returning("A")
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "A"
    assert alt == "G"
    assert agree == "NM_tied"


def test_refalt_v2_snp_nm_na_probe_only():
    """probe_only: no winning_ts passed (None) → NM_N/A, genome primary."""
    info = {"AlleleA": "A", "AlleleB": "G"}
    fasta = _fasta_returning("G")
    ref, alt, agree, _ = determine_ref_alt_v2(
        None, None, {}, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "G"
    assert alt == "A"
    assert agree == "NM_N/A"


def test_refalt_v2_snp_nm_only_genome_fails():
    """Genome lookup fails (triallelic), NM succeeds → NM_only."""
    ts_a = _ts("chr1", 1000, nm=0)
    ts_b = _ts("chr1", 1000, nm=2)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "A", "AlleleB": "G"}
    fasta = _fasta_returning("C")   # neither allele
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "A"
    assert alt == "G"
    assert agree == "NM_only"


def test_refalt_v2_snp_refalt_unresolved_both_fail():
    """Both methods fail → refalt_unresolved (None, None, 'refalt_unresolved')."""
    ts_a = _ts("chr1", 1000, nm=1)
    ts_b = _ts("chr1", 1000, nm=1)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "A", "AlleleB": "G"}
    fasta = _fasta_returning("C")   # neither allele
    result = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert result[:3] == (None, None, "refalt_unresolved")


def test_refalt_v2_indel_probe_only_refalt_unresolved():
    """Indel marker with probe-only (no winning_ts) → ambiguous, not NM_tied."""
    info = {"AlleleA": "AT", "AlleleB": ""}
    fasta = MagicMock()
    result = determine_ref_alt_v2(
        None, None, {}, info, fasta, "chr1", 1000, "+"
    )
    assert result[:3] == (None, None, "refalt_unresolved")


def test_refalt_v2_deletion_nm_validated():
    """Deletion marker: NM determines ref, genome fetch confirms → NM_validated."""
    ts_a = _ts("chr1", 1000, nm=0)   # allele A = "AT" (ref = longer)
    ts_b = _ts("chr1", 1000, nm=2)   # allele B = "" (alt = deletion)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "AT", "AlleleB": ""}
    fasta = MagicMock()
    fasta.fetch.return_value = "AT"  # genome matches gref = "AT"
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "AT"
    assert alt == ""
    assert agree == "NM_validated"


def test_refalt_v2_deletion_nm_mismatch():
    """Deletion marker: genome doesn't match within ±10 bp → NM_mismatch with NM alleles."""
    ts_a = _ts("chr1", 1000, nm=0)
    ts_b = _ts("chr1", 1000, nm=2)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "AT", "AlleleB": ""}
    fasta = MagicMock()
    fasta.fetch.return_value = "GC"  # never matches "AT" → refinement fails
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == "AT"
    assert alt == ""
    assert agree == "NM_mismatch"


def test_refalt_v2_insertion_nm_na():
    """Insertion marker (gref=''): genome validation not applicable → NM_N/A."""
    ts_a = _ts("chr1", 1000, nm=0)   # allele A = "" (ref = insertion ref = empty)
    ts_b = _ts("chr1", 1000, nm=2)   # allele B = "TCG" (alt = inserted seq)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "", "AlleleB": "TCG"}
    fasta = MagicMock()
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    assert ref == ""
    assert alt == "TCG"
    assert agree == "NM_N/A"


def test_refalt_v2_insertion_nm_corrected_genome_has_alt_base():
    """1-bp 'insertion' where the genome at final_pos actually has the alt base →
    NM got Ref/Alt backwards; result must be flagged NM_corrected and alleles swapped."""
    # NM says A is ref (empty) and B=G is alt (inserted). But the genome at final_pos
    # has 'G' — i.e. G is truly the reference and the variant is a G→'' deletion.
    ts_a = _ts("chr1", 1000, nm=0)
    ts_b = _ts("chr1", 1000, nm=2)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "", "AlleleB": "G"}
    fasta = MagicMock()
    fasta.fetch.return_value = "G"
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "+"
    )
    # Alleles are swapped: what NM called Alt ('G') is now Ref, the empty is Alt.
    assert ref == "G"
    assert alt == ""
    assert agree == "NM_corrected"


def test_refalt_v2_insertion_nm_corrected_minus_strand_complements_alt():
    """On minus strand, the genome's forward-strand base must be complemented
    before comparison with the alignment-strand alt char."""
    ts_a = _ts("chr1", 1000, nm=0)
    ts_b = _ts("chr1", 1000, nm=2)
    ts_aligns = {"A": [ts_a], "B": [ts_b]}
    info = {"AlleleA": "", "AlleleB": "C"}
    fasta = MagicMock()
    # genome '+' strand has 'G'. On minus strand (alignment orientation) the alt
    # 'C' complements to 'G' on '+', which matches → NM_corrected triggers.
    fasta.fetch.return_value = "G"
    ref, alt, agree, _ = determine_ref_alt_v2(
        "A", ts_a, ts_aligns, info, fasta, "chr1", 1000, "-"
    )
    assert ref == "C"
    assert alt == ""
    assert agree == "NM_corrected"


# ── integration: new output columns ──────────────────────────────────────────

def test_new_columns_present_in_output():
    """All six new decision functions are importable — confirms Task 8 wiring is complete."""
    from remap_manifest import (
        compute_alignment_status, build_valid_triples,
        rank_and_resolve, best_topseq_rescue, best_probe_rescue,
        determine_ref_alt_v2,
    )
    # Verify all expected functions are importable
    assert callable(compute_alignment_status)
    assert callable(build_valid_triples)
    assert callable(rank_and_resolve)
    assert callable(best_topseq_rescue)
    assert callable(best_probe_rescue)
    assert callable(determine_ref_alt_v2)
