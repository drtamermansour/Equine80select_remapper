import os
import pytest
import pandas as pd
from benchmark_compare import (
    normalise_chr, classify_marker, write_diff,
    stratify_by_coord_delta,
    parse_snp_alleles, parse_topseq_alleles,
    alleles_match_snp, check_flanking_context,
    classify_explanatory,
    compute_qc_impact, QC_STAGE_ORDER,
)


# ── normalise_chr ─────────────────────────────────────────────────────────────

def test_normalise_chr_x_alias():
    assert normalise_chr("X_NC_009175.3") == "X"

def test_normalise_chr_x_unchanged():
    assert normalise_chr("X") == "X"

def test_normalise_chr_autosome_unchanged():
    assert normalise_chr("1") == "1"
    assert normalise_chr("31") == "31"

def test_normalise_chr_y_unchanged():
    assert normalise_chr("Y") == "Y"

def test_normalise_chr_zero_unchanged():
    assert normalise_chr("0") == "0"


# ── classify_marker ───────────────────────────────────────────────────────────

def _mrow(chr="1", pos=1000, strand="+"):
    return {"manifest_chr": chr, "manifest_pos": pos, "manifest_strand": strand}

def _rrow(chr="1", pos=1000, strand="+", status="mapped"):
    return {
        "remapped_chr": chr, "remapped_pos": pos,
        "remapped_strand": strand, "remapped_status": status,
    }

def test_classify_correct():
    assert classify_marker(_mrow(), _rrow()) == "correct"

def test_classify_correct_x_alias():
    m = _mrow(chr="X_NC_009175.3")
    r = _rrow(chr="X")
    assert classify_marker(m, r) == "correct"

def test_classify_coord_correct_strand_wrong():
    assert classify_marker(_mrow(), _rrow(strand="-")) == "coord_correct_strand_wrong"

def test_classify_coord_off():
    assert classify_marker(_mrow(pos=1000), _rrow(pos=1051)) == "coord_off"

def test_classify_wrong_chr():
    assert classify_marker(_mrow(chr="1"), _rrow(chr="2")) == "wrong_chr"

def test_classify_unmapped_chr0():
    assert classify_marker(_mrow(), _rrow(chr="0")) == "unmapped"

def test_classify_unmapped_strand_na():
    assert classify_marker(_mrow(), _rrow(strand="N/A")) == "unmapped"

def test_classify_ambiguous():
    assert classify_marker(_mrow(), _rrow(status="ambiguous")) == "ambiguous"

def test_classify_unmapped_takes_priority_over_ambiguous():
    assert classify_marker(_mrow(), _rrow(chr="0", status="ambiguous")) == "unmapped"

def test_classify_coord_off_non_numeric_remapped_pos():
    # Non-numeric remapped position → coord_off (caught by ValueError)
    assert classify_marker(_mrow(pos=1000), _rrow(pos="invalid")) == "coord_off"

def test_classify_coord_off_none_remapped_pos():
    # None remapped position → coord_off
    assert classify_marker(_mrow(pos=1000), _rrow(pos=None)) == "coord_off"


def test_classify_uses_probe_strand_when_present():
    """When remapped_probe_strand is '+' or '-', it takes precedence over TopSeq strand."""
    m = _mrow(strand="+")
    # TopSeq strand says - but probe strand says + → probe strand wins → correct
    r = _rrow(strand="-")
    r["remapped_probe_strand"] = "+"
    assert classify_marker(m, r) == "correct"


def test_classify_probe_strand_mismatch_is_strand_wrong():
    """Probe strand differs from manifest strand (RefStrand) → coord_correct_strand_wrong."""
    m = _mrow(strand="+")
    r = _rrow(strand="+")  # TopSeq strand matches
    r["remapped_probe_strand"] = "-"  # but probe strand disagrees
    assert classify_marker(m, r) == "coord_correct_strand_wrong"


def test_classify_probe_strand_na_exempts_topseq_only():
    """remapped_probe_strand == 'N/A' → strand check skipped (topseq_only markers)."""
    m = _mrow(strand="+")
    r = _rrow(strand="-")
    r["remapped_probe_strand"] = "N/A"
    assert classify_marker(m, r) == "correct"


def test_classify_falls_back_to_topseq_strand_when_probe_strand_absent():
    """Backward compat: missing remapped_probe_strand key → use remapped_strand."""
    m = _mrow(strand="+")
    r = _rrow(strand="+")    # no probe_strand key
    assert classify_marker(m, r) == "correct"
    m2 = _mrow(strand="+")
    r2 = _rrow(strand="-")
    assert classify_marker(m2, r2) == "coord_correct_strand_wrong"


# ── parse_snp_alleles ─────────────────────────────────────────────────────────

def test_parse_snp_alleles_standard_snp():
    assert parse_snp_alleles("[A/G]") == ("A", "G")

def test_parse_snp_alleles_indel_returns_none():
    """[D/I] is not a real sequence pair — return None."""
    assert parse_snp_alleles("[D/I]") is None
    assert parse_snp_alleles("[I/D]") is None

def test_parse_snp_alleles_malformed_returns_none():
    assert parse_snp_alleles("no_brackets") is None
    assert parse_snp_alleles("") is None
    assert parse_snp_alleles(None) is None


# ── parse_topseq_alleles ──────────────────────────────────────────────────────

def test_parse_topseq_alleles_snp():
    """TopGenomicSeq with SNP: returns (prefix, A, B, suffix)."""
    pre, a, b, post = parse_topseq_alleles("ACGTACGT[A/G]TACGTACG")
    assert (pre, a, b, post) == ("ACGTACGT", "A", "G", "TACGTACG")

def test_parse_topseq_alleles_deletion():
    """Deletion — 'B' is empty string."""
    pre, a, b, post = parse_topseq_alleles("ACGTACGT[CCCGGG/-]TACGTACG")
    assert (pre, a, b, post) == ("ACGTACGT", "CCCGGG", "", "TACGTACG")

def test_parse_topseq_alleles_insertion():
    """Insertion — 'A' is empty string."""
    pre, a, b, post = parse_topseq_alleles("ACGTACGT[-/CCCGGG]TACGTACG")
    assert (pre, a, b, post) == ("ACGTACGT", "", "CCCGGG", "TACGTACG")

def test_parse_topseq_alleles_malformed_returns_none():
    assert parse_topseq_alleles("no_brackets") is None
    assert parse_topseq_alleles("") is None


# ── alleles_match_snp ─────────────────────────────────────────────────────────

def test_alleles_match_snp_direct_on_plus_strand():
    """Remapped + strand, alleles match manifest set → True."""
    # Remapped Ref=A, Alt=G on + strand; manifest SNP=[A/G]
    assert alleles_match_snp("A", "G", "+", ("A", "G")) is True

def test_alleles_match_snp_swap_ok():
    """Order of manifest and remapped alleles doesn't matter (set comparison)."""
    assert alleles_match_snp("G", "A", "+", ("A", "G")) is True

def test_alleles_match_snp_minus_strand_needs_rc():
    """Remapped on - strand: Ref/Alt must RC to match manifest."""
    # Remapped Ref=T, Alt=C on - strand → fwd_Ref=A, fwd_Alt=G; matches {A,G}
    assert alleles_match_snp("T", "C", "-", ("A", "G")) is True

def test_alleles_match_snp_ambiguous_pair_at_passes():
    """A/T SNP: {A,T} equals its own complement, so either strand matches."""
    assert alleles_match_snp("A", "T", "+", ("A", "T")) is True
    assert alleles_match_snp("A", "T", "-", ("A", "T")) is True   # RC = {T,A}

def test_alleles_match_snp_mismatch():
    """Different alleles that do not match directly *or* via complement → False."""
    # Remapped {A, G} → direct set {A, G}, complement set {T, C}
    # Manifest {A, C} is not equal to either → False
    assert alleles_match_snp("A", "G", "+", ("A", "C")) is False

def test_alleles_match_snp_complement_only_match():
    """Manifest and remapped alleles only match after both-complement — True.

    Covers the case where RefStrand says + but the SNP column is in the opposite
    convention: {complement(Ref), complement(Alt)} == {SNP_A, SNP_B}.
    """
    # Remapped Ref=A, Alt=G (stored on + strand); manifest SNP=[T/C].
    # {A, G} != {T, C}; but {complement(A), complement(G)} == {T, C} → True.
    assert alleles_match_snp("A", "G", "+", ("T", "C")) is True


# ── check_flanking_context ────────────────────────────────────────────────────

class _FakeFasta:
    """Minimal stub of pysam.FastaFile for tests."""
    def __init__(self, seqs):
        self._seqs = seqs     # {chrom: full_sequence}
    def fetch(self, chrom, start, end):
        return self._seqs[chrom][start:end]


def test_check_flanking_context_forward_match():
    """Genome left == PREFIX suffix, genome right == SUFFIX prefix (forward)."""
    # Place the variant at 1-based position 21; PREFIX=10bp, SUFFIX=10bp
    seq = "ACGTACGTAC" + "N" + "TACGTACGTA"   # 21 bp, variant is 'N' at index 10
    fasta = _FakeFasta({"chr1": seq})
    fwd, rev = check_flanking_context(fasta, "chr1", mapinfo=11, prefix="ACGTACGTAC",
                                       suffix="TACGTACGTA", allele_len=1, flank_len=10)
    assert fwd is True
    assert rev is False


def test_check_flanking_context_reverse_match():
    """Genome flanking matches RC(SUFFIX) and RC(PREFIX) — reverse orientation."""
    import importlib
    from remap_manifest import reverse_complement
    prefix = "ACGTACGTAC"
    suffix = "TACGTACGTA"
    # Build genome that contains RC(SUFFIX) before, RC(PREFIX) after the variant
    seq = reverse_complement(suffix) + "N" + reverse_complement(prefix)
    fasta = _FakeFasta({"chr1": seq})
    fwd, rev = check_flanking_context(fasta, "chr1", mapinfo=11, prefix=prefix,
                                       suffix=suffix, allele_len=1, flank_len=10)
    assert fwd is False
    assert rev is True


def test_check_flanking_context_no_match():
    """Neither orientation matches — both False."""
    fasta = _FakeFasta({"chr1": "G" * 100})
    fwd, rev = check_flanking_context(fasta, "chr1", mapinfo=50, prefix="ACGTACGTAC",
                                       suffix="TACGTACGTA", allele_len=1, flank_len=10)
    assert fwd is False
    assert rev is False


# ── classify_explanatory ──────────────────────────────────────────────────────

def test_classify_explanatory_manifest_strand_wrong():
    """context_forward True + remapped_strand != RefStrand → verdict manifest_strand_wrong."""
    row = {
        "context_forward": True, "context_reverse": False,
        "remapped_strand": "+", "manifest_strand": "-",
        "coord_ok": True, "is_ambiguous_snp": False, "is_probe_only": False,
        "is_indel": False,
        "deletion_seq_ok": None, "insertion_absent": None,
    }
    assert classify_explanatory(row) == "manifest_strand_wrong"


def test_classify_explanatory_pipeline_wrong_locus():
    """Context fails at remapped MapInfo → pipeline_wrong_locus."""
    row = {
        "context_forward": False, "context_reverse": False,
        "remapped_strand": "+", "manifest_strand": "+",
        "coord_ok": True, "is_ambiguous_snp": False, "is_probe_only": False,
        "is_indel": False,
        "deletion_seq_ok": None, "insertion_absent": None,
    }
    assert classify_explanatory(row) == "pipeline_wrong_locus"


def test_classify_explanatory_pipeline_wrong_strand():
    """context_reverse True, context_forward False, strand disagrees → pipeline_wrong_strand."""
    row = {
        "context_forward": False, "context_reverse": True,
        "remapped_strand": "+", "manifest_strand": "-",
        "coord_ok": True, "is_ambiguous_snp": False, "is_probe_only": False,
        "is_indel": False,
        "deletion_seq_ok": None, "insertion_absent": None,
    }
    assert classify_explanatory(row) == "pipeline_wrong_strand"


def test_classify_explanatory_ambiguous_snp():
    row = {
        "context_forward": True, "context_reverse": True,
        "remapped_strand": "+", "manifest_strand": "+",
        "coord_ok": True, "is_ambiguous_snp": True, "is_probe_only": False,
        "is_indel": False,
        "deletion_seq_ok": None, "insertion_absent": None,
    }
    assert classify_explanatory(row) == "ambiguous_snp"


def test_classify_explanatory_probe_only_inconclusive():
    row = {
        "context_forward": False, "context_reverse": False,
        "remapped_strand": "+", "manifest_strand": "+",
        "coord_ok": True, "is_ambiguous_snp": False, "is_probe_only": True,
        "is_indel": False,
        "deletion_seq_ok": None, "insertion_absent": None,
    }
    assert classify_explanatory(row) == "probe_only_inconclusive"


def test_classify_explanatory_unresolved():
    row = {
        "context_forward": True, "context_reverse": False,
        "remapped_strand": "+", "manifest_strand": "+",
        "coord_ok": False, "is_ambiguous_snp": False, "is_probe_only": False,
        "is_indel": False,
        "deletion_seq_ok": None, "insertion_absent": None,
    }
    # context matches at our coord, our coord != manifest coord → manifest_coord_wrong
    assert classify_explanatory(row) == "manifest_coord_wrong"


# ── compute_qc_impact ─────────────────────────────────────────────────────────

def _impact_fixture():
    """Small synthetic result_df for QC impact tests.

    Rows:
      - 3 passed all QC + correct
      - 1 passed all QC + coord_off  (false negative)
      - 2 removed at stage_5 + correct (false positive x2)
      - 1 removed at stage_5 + coord_off
      - 1 removed at stage_10 + unmapped
    Total 8 rows.
    """
    return pd.DataFrame({
        "Name":         [f"m{i}" for i in range(8)],
        "why_filtered": ["", "", "", "",
                         "stage_5_refalt_conf", "stage_5_refalt_conf",
                         "stage_5_refalt_conf", "stage_10_polymorphic"],
        "result":       ["correct", "correct", "correct", "coord_off",
                         "correct", "correct",
                         "coord_off", "unmapped"],
    })


def test_compute_qc_impact_passed_counts():
    impact = compute_qc_impact(_impact_fixture())
    assert impact["passed_n"] == 4
    assert impact["passed_correct"] == 3
    assert impact["passed_accuracy_pct"] == pytest.approx(75.0)


def test_compute_qc_impact_per_stage_precision():
    impact = compute_qc_impact(_impact_fixture())
    per_stage = dict((s, (n, nc, pct)) for s, n, nc, pct in impact["per_stage"])
    # stage_5: 3 removed, 1 non-correct (the coord_off one) → precision 33.3%
    assert per_stage["stage_5_refalt_conf"] == pytest.approx((3, 1, 33.333333), abs=0.01)
    # stage_10: 1 removed, 1 non-correct → precision 100%
    assert per_stage["stage_10_polymorphic"] == pytest.approx((1, 1, 100.0), abs=0.01)


def test_compute_qc_impact_fp_and_fn_identification():
    impact = compute_qc_impact(_impact_fixture())
    # FP: 2 correct markers removed at stage_5
    assert len(impact["fp_df"]) == 2
    assert set(impact["fp_df"]["Name"]) == {"m4", "m5"}
    # FN: 1 non-correct marker surviving (coord_off at m3)
    assert len(impact["fn_df"]) == 1
    assert list(impact["fn_df"]["Name"]) == ["m3"]


def test_compute_qc_impact_cumulative_accuracy_monotonic_to_final():
    """Cumulative passing-set ends at the 'all QC applied' count."""
    impact = compute_qc_impact(_impact_fixture())
    final_label, final_n, final_c, final_acc = impact["cumulative"][-1]
    assert final_label == QC_STAGE_ORDER[-1]   # stage_11_ambiguous_snp (no rejections in fixture)
    assert final_n == impact["passed_n"]
    assert final_c == impact["passed_correct"]


def test_compute_qc_impact_cumulative_accuracy_rises_when_stage_removes_only_errors():
    """When a stage removes only non-correct markers, accuracy goes up."""
    df = pd.DataFrame({
        "Name": ["a", "b", "c", "d"],
        "why_filtered": ["stage_10_polymorphic", "stage_10_polymorphic", "", ""],
        "result": ["coord_off", "unmapped", "correct", "correct"],
    })
    impact = compute_qc_impact(df)
    # Before QC: 4 markers, 2 correct → 50%
    assert impact["cumulative"][0] == pytest.approx(("(before QC)", 4, 2, 50.0), abs=0.01)
    # After stage_10 removes the 2 non-correct ones: 2 markers, 2 correct → 100%
    after_s10 = next(t for t in impact["cumulative"] if t[0] == "stage_10_polymorphic")
    assert after_s10 == pytest.approx(("stage_10_polymorphic", 2, 2, 100.0), abs=0.01)


# ── load_manifest ─────────────────────────────────────────────────────────────

import textwrap, io
from benchmark_compare import load_manifest


MANIFEST_CSV = textwrap.dedent("""\
    Illumina, Inc.,,,,,,,,,,,,,,,,,,,
    [Heading],,,,,,,,,,,,,,,,,,,,
    Descriptor File Name,test.bpm,,,,,,,,,,,,,,,,,,,
    Assay Format,Infinium HTS,,,,,,,,,,,,,,,,,,,
    Date Manufactured,1/1/2025,,,,,,,,,,,,,,,,,,,
    Loci Count ,6,,,,,,,,,,,,,,,,,,,
    [Assay],,,,,,,,,,,,,,,,,,,,
    IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand
    id1,SNP_auto,TOP,[A/G],,,,,3,1,1000,diploid,Equus caballus,Equcab,3,TOP,seq,seq,1,3,+
    id2,SNP_x,TOP,[A/G],,,,,3,X,2000,diploid,Equus caballus,Equcab,3,TOP,seq,seq,1,3,+
    id3,SNP_x_alias,TOP,[A/G],,,,,3,X_NC_009175.3,3000,diploid,Equus caballus,Equcab,3,TOP,seq,seq,1,3,-
    id4,SNP_y,TOP,[A/G],,,,,3,Y,4000,diploid,Equus caballus,Equcab,3,TOP,seq,seq,1,3,+
    id5,SNP_unplaced,TOP,[A/G],,,,,3,0,0,diploid,Equus caballus,Equcab,3,TOP,seq,seq,1,3,+
    id6,SNP_auto2,TOP,[A/G],,,,,3,2,5000,diploid,Equus caballus,Equcab,3,TOP,seq,seq,1,3,-
    [Controls],,,,,,,,,,,,,,,,,,,,
    control1,ctrl,,,,,,,,,,,,,,,,,,,
""")


@pytest.fixture
def manifest_file(tmp_path):
    p = tmp_path / "test_manifest.csv"
    p.write_text(MANIFEST_CSV)
    return str(p)


def test_load_manifest_main_scope(manifest_file):
    main_df, chry_df, chr0_df = load_manifest(manifest_file)
    assert set(main_df["Name"]) == {"SNP_auto", "SNP_x", "SNP_x_alias", "SNP_auto2"}

def test_load_manifest_chry(manifest_file):
    main_df, chry_df, chr0_df = load_manifest(manifest_file)
    assert list(chry_df["Name"]) == ["SNP_y"]

def test_load_manifest_chr0(manifest_file):
    main_df, chry_df, chr0_df = load_manifest(manifest_file)
    assert list(chr0_df["Name"]) == ["SNP_unplaced"]

def test_load_manifest_x_alias_normalised(manifest_file):
    main_df, chry_df, chr0_df = load_manifest(manifest_file)
    x_rows = main_df[main_df["Name"] == "SNP_x_alias"]
    assert x_rows.iloc[0]["manifest_chr"] == "X"

def test_load_manifest_columns(manifest_file):
    main_df, _, _ = load_manifest(manifest_file)
    assert set(main_df.columns) >= {"Name", "manifest_chr", "manifest_pos", "manifest_strand"}

def test_load_manifest_controls_excluded(manifest_file):
    main_df, chry_df, chr0_df = load_manifest(manifest_file)
    all_names = set(main_df["Name"]) | set(chry_df["Name"]) | set(chr0_df["Name"])
    assert "control1" not in all_names


# ── load_remapped ────────────────────────────────────────────────────────────

from benchmark_compare import load_remapped, compare_all


REMAPPED_CSV = textwrap.dedent("""\
    Name,Chr_equCab3,MapInfo_equCab3,Strand_equCab3,MappingStatus_equCab3
    SNP_auto,1,1000,+,mapped
    SNP_x,X,2000,+,mapped
    SNP_x_alias,X,3000,-,mapped
    SNP_auto2,3,5000,-,mapped
    SNP_wrong_chr,5,9999,+,mapped
    SNP_coord_off,1,1100,+,mapped
    SNP_strand_wrong,1,7000,+,mapped
    SNP_unmapped,0,0,N/A,unmapped
    SNP_ambiguous,1,8000,+,ambiguous
""")

MANIFEST_MAIN = textwrap.dedent("""\
    Name,manifest_chr,manifest_pos,manifest_strand
    SNP_auto,1,1000,+
    SNP_x,X,2000,+
    SNP_x_alias,X,3000,-
    SNP_auto2,2,5000,-
    SNP_wrong_chr,1,9999,+
    SNP_coord_off,1,1000,+
    SNP_strand_wrong,1,7000,-
    SNP_unmapped,1,6000,+
    SNP_ambiguous,1,8000,+
    SNP_missing_from_remapped,1,9000,+
""")


@pytest.fixture
def remapped_file(tmp_path):
    p = tmp_path / "remapped.csv"
    p.write_text(REMAPPED_CSV)
    return str(p)


@pytest.fixture
def manifest_main_df():
    return pd.read_csv(io.StringIO(MANIFEST_MAIN), dtype={"manifest_chr": str, "manifest_pos": int})


def test_load_remapped_columns(remapped_file):
    df = load_remapped(remapped_file, "equCab3")
    assert set(df.columns) >= {"Name", "remapped_chr", "remapped_pos",
                                "remapped_strand", "remapped_status"}

def test_load_remapped_x_normalised(remapped_file):
    df = load_remapped(remapped_file, "equCab3")
    assert "X" in df["remapped_chr"].values

def test_compare_all_correct(manifest_main_df, remapped_file):
    remapped_df = load_remapped(remapped_file, "equCab3")
    result = compare_all(manifest_main_df, remapped_df)
    assert result.loc[result["Name"] == "SNP_auto", "result"].iloc[0] == "correct"

def test_compare_all_wrong_chr(manifest_main_df, remapped_file):
    remapped_df = load_remapped(remapped_file, "equCab3")
    result = compare_all(manifest_main_df, remapped_df)
    assert result.loc[result["Name"] == "SNP_auto2", "result"].iloc[0] == "wrong_chr"

def test_compare_all_coord_off(manifest_main_df, remapped_file):
    remapped_df = load_remapped(remapped_file, "equCab3")
    result = compare_all(manifest_main_df, remapped_df)
    row = result.loc[result["Name"] == "SNP_coord_off"].iloc[0]
    assert row["result"] == "coord_off"
    assert row["coord_offset"] == 100

def test_compare_all_strand_wrong(manifest_main_df, remapped_file):
    remapped_df = load_remapped(remapped_file, "equCab3")
    result = compare_all(manifest_main_df, remapped_df)
    assert result.loc[result["Name"] == "SNP_strand_wrong", "result"].iloc[0] == "coord_correct_strand_wrong"

def test_compare_all_missing_from_remapped_is_unmapped(manifest_main_df, remapped_file):
    remapped_df = load_remapped(remapped_file, "equCab3")
    result = compare_all(manifest_main_df, remapped_df)
    assert result.loc[result["Name"] == "SNP_missing_from_remapped", "result"].iloc[0] == "unmapped"

def test_compare_all_coord_offset_blank_for_wrong_chr(manifest_main_df, remapped_file):
    remapped_df = load_remapped(remapped_file, "equCab3")
    result = compare_all(manifest_main_df, remapped_df)
    offset = result.loc[result["Name"] == "SNP_wrong_chr", "coord_offset"].iloc[0]
    assert pd.isna(offset)


from benchmark_compare import write_tsv, write_report, load_baseline


def _make_result_df(rows):
    """Helper: build a minimal result DataFrame from a list of dicts."""
    return pd.DataFrame(rows)


def test_write_tsv_creates_file(tmp_path):
    df = _make_result_df([
        {"Name": "SNP1", "manifest_chr": "1", "manifest_pos": 1000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 1000,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "correct", "coord_offset": None},
    ])
    out = str(tmp_path / "out.tsv")
    write_tsv(df, out)
    assert os.path.exists(out)


def test_write_tsv_columns(tmp_path):
    df = _make_result_df([
        {"Name": "SNP1", "manifest_chr": "1", "manifest_pos": 1000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 1000,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "correct", "coord_offset": None},
    ])
    out = str(tmp_path / "out.tsv")
    write_tsv(df, out)
    result = pd.read_csv(out, sep="\t")
    expected_cols = {"Name", "manifest_chr", "manifest_pos", "manifest_strand",
                     "remapped_chr", "remapped_pos", "remapped_strand",
                     "remapped_status", "result", "coord_offset"}
    assert expected_cols.issubset(set(result.columns))


def test_write_report_contains_headline(tmp_path):
    df = _make_result_df([
        {"Name": "SNP1", "manifest_chr": "1", "manifest_pos": 1000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 1000,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "correct", "coord_offset": None},
        {"Name": "SNP2", "manifest_chr": "1", "manifest_pos": 2000,
         "manifest_strand": "+", "remapped_chr": "0", "remapped_pos": 0,
         "remapped_strand": "N/A", "remapped_status": "unmapped",
         "result": "unmapped", "coord_offset": None},
    ])
    chry_df = _make_result_df([])
    chr0_df = _make_result_df([])
    out = str(tmp_path / "report.txt")
    write_report(df, chry_df, chr0_df, out, assembly="equCab3",
                 manifest_path="test.csv", remapped_path="test_remapped.csv")
    text = open(out).read()
    assert "HEADLINE COUNTS" in text
    assert "correct" in text
    assert "unmapped" in text


def test_load_baseline_roundtrip(tmp_path):
    df = _make_result_df([
        {"Name": "SNP1", "manifest_chr": "1", "manifest_pos": 1000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 1000,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "correct", "coord_offset": None},
    ])
    out = str(tmp_path / "baseline.tsv")
    write_tsv(df, out)
    loaded = load_baseline(out)
    assert list(loaded["Name"]) == ["SNP1"]
    assert list(loaded["result"]) == ["correct"]


# ── stratify_by_coord_delta ───────────────────────────────────────────────────

def _result_df_with_delta(rows):
    """Build a minimal result DataFrame that includes coord_delta."""
    return pd.DataFrame(rows)


def test_stratify_returns_none_when_column_absent():
    df = pd.DataFrame([{"Name": "S1", "result": "correct"}])
    assert stratify_by_coord_delta(df) is None


def test_stratify_basic_accuracy():
    df = _result_df_with_delta([
        {"Name": "S1", "result": "correct",                    "coord_delta": 0},
        {"Name": "S2", "result": "correct",                    "coord_delta": 0},
        {"Name": "S3", "result": "coord_correct_strand_wrong", "coord_delta": 0},
        {"Name": "S4", "result": "coord_off",                  "coord_delta": 1},
        {"Name": "S5", "result": "wrong_chr",                  "coord_delta": 5},
        {"Name": "S6", "result": "unmapped",                   "coord_delta": -1},
    ])
    rows = stratify_by_coord_delta(df)
    assert rows is not None

    by_label = {r["label"]: r for r in rows}

    d0 = by_label["delta = 0"]
    assert d0["n"] == 3
    assert d0["coord_accurate"] == 3   # correct + coord_correct_strand_wrong
    assert d0["correct"] == 2

    d1 = by_label["delta = 1"]
    assert d1["n"] == 1
    assert d1["coord_accurate"] == 0

    d210 = by_label["delta = 2-10"]
    assert d210["n"] == 1
    assert d210["coord_accurate"] == 0

    dm1 = by_label["delta = -1"]
    assert dm1["n"] == 1
    assert dm1["coord_accurate"] == 0


def test_stratify_empty_bucket_excluded():
    # No delta > 10 markers — that bucket should not appear
    df = _result_df_with_delta([
        {"Name": "S1", "result": "correct", "coord_delta": 0},
    ])
    rows = stratify_by_coord_delta(df)
    labels = [r["label"] for r in rows]
    assert "delta > 10" not in labels


def test_stratify_delta_gt10_bucket():
    df = _result_df_with_delta([
        {"Name": "S1", "result": "correct", "coord_delta": 0},
        {"Name": "S2", "result": "coord_off", "coord_delta": 50},
    ])
    rows = stratify_by_coord_delta(df)
    by_label = {r["label"]: r for r in rows}
    assert "delta > 10" in by_label
    assert by_label["delta > 10"]["n"] == 1
    assert by_label["delta > 10"]["coord_accurate"] == 0


# ── load_remapped with coord_delta ────────────────────────────────────────────

REMAPPED_CSV_WITH_DELTA = textwrap.dedent("""\
    Name,Chr_equCab3,MapInfo_equCab3,Strand_equCab3,MappingStatus_equCab3,CoordDelta_equCab3
    SNP_auto,1,1000,+,mapped,0
    SNP_coord_off,1,1100,+,mapped,1
    SNP_unmapped,0,0,N/A,unmapped,-1
""")


@pytest.fixture
def remapped_file_with_delta(tmp_path):
    p = tmp_path / "remapped_delta.csv"
    p.write_text(REMAPPED_CSV_WITH_DELTA)
    return str(p)


def test_load_remapped_loads_coord_delta_when_present(remapped_file_with_delta):
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_with_delta, "equCab3")
    assert "coord_delta" in df.columns
    row = df.loc[df["Name"] == "SNP_auto"].iloc[0]
    assert row["coord_delta"] == 0


def test_load_remapped_no_coord_delta_column_when_absent(remapped_file):
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file, "equCab3")
    assert "coord_delta" not in df.columns


# ── write_report includes coord_delta section when column present ─────────────

def test_write_report_includes_coord_delta_section(tmp_path):
    df = _make_result_df([
        {"Name": "S1", "manifest_chr": "1", "manifest_pos": 1000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 1000,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "correct", "coord_offset": None, "coord_delta": 0},
        {"Name": "S2", "manifest_chr": "1", "manifest_pos": 2000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 2001,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "coord_off", "coord_offset": 1, "coord_delta": 1},
    ])
    chry_df = _make_result_df([])
    chr0_df = _make_result_df([])
    out = str(tmp_path / "report.txt")
    write_report(df, chry_df, chr0_df, out, assembly="equCab3",
                 manifest_path="test.csv", remapped_path="test_remapped.csv")
    text = open(out).read()
    assert "ACCURACY STRATIFIED BY COORD_DELTA" in text
    assert "delta = 0" in text
    assert "delta = 1" in text


def test_write_report_omits_coord_delta_section_when_column_absent(tmp_path):
    df = _make_result_df([
        {"Name": "S1", "manifest_chr": "1", "manifest_pos": 1000,
         "manifest_strand": "+", "remapped_chr": "1", "remapped_pos": 1000,
         "remapped_strand": "+", "remapped_status": "mapped",
         "result": "correct", "coord_offset": None},
    ])
    chry_df = _make_result_df([])
    chr0_df = _make_result_df([])
    out = str(tmp_path / "report.txt")
    write_report(df, chry_df, chr0_df, out, assembly="equCab3",
                 manifest_path="test.csv", remapped_path="test_remapped.csv")
    text = open(out).read()
    assert "ACCURACY STRATIFIED BY COORD_DELTA" not in text


import subprocess


def test_integration_runs_without_error(tmp_path, results_dir):
    """Smoke test: run the script against real data and check output files exist."""
    result = subprocess.run(
        [
            "python", "scripts/benchmark_compare.py",
            "--manifest",   "backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv",
            "--remapped",   os.path.join(results_dir, "remapping",
                            "Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv"),
            "--assembly",   "equCab3",
            "--output-dir", str(tmp_path),
        ],
        cwd="/home/tahmed/Equine80select_remapper",
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    files = os.listdir(tmp_path)
    assert any(f.endswith(".tsv") and "chrY" not in f and "chr0" not in f for f in files), \
        f"Main TSV missing. Files: {files}"
    assert any("chrY" in f for f in files),  f"chrY TSV missing. Files: {files}"
    assert any("chr0" in f for f in files),  f"chr0 TSV missing. Files: {files}"
    assert any("report" in f for f in files), f"Report missing. Files: {files}"


def test_integration_correct_count(tmp_path, results_dir):
    """Correct marker count should be well above 80,000 for the current run."""
    subprocess.run(
        [
            "python", "scripts/benchmark_compare.py",
            "--manifest",   "backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv",
            "--remapped",   os.path.join(results_dir, "remapping",
                            "Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv"),
            "--assembly",   "equCab3",
            "--output-dir", str(tmp_path),
        ],
        cwd="/home/tahmed/Equine80select_remapper",
        capture_output=True,
    )
    tsv_files = [f for f in os.listdir(tmp_path)
                 if f.endswith(".tsv") and "chrY" not in f and "chr0" not in f]
    df = pd.read_csv(os.path.join(tmp_path, tsv_files[0]), sep="\t")
    correct_count = (df["result"] == "correct").sum()
    assert correct_count > 80_000, f"Expected >80k correct, got {correct_count}"


# ── write_diff ────────────────────────────────────────────────────────────────

def test_write_diff_creates_file(tmp_path):
    """write_diff writes a file listing category transitions."""
    curr = pd.DataFrame({
        "Name": ["A", "B", "C"],
        "result": ["correct", "unmapped", "coord_off"],
    })
    base = pd.DataFrame({
        "Name": ["A", "B", "C"],
        "result": ["unmapped", "unmapped", "correct"],
    })
    diff_path = str(tmp_path / "diff.txt")
    write_diff(curr, base, diff_path)
    assert os.path.exists(diff_path)
    content = open(diff_path).read()
    assert "Markers that changed category: 2" in content
    assert "unmapped" in content
    assert "correct" in content


def test_integration_baseline_produces_diff_file(tmp_path, results_dir):
    """When --baseline is provided, a _diff.txt file is created alongside the report."""
    run1 = tmp_path / "run1"
    run1.mkdir()
    remapped_csv = os.path.join(results_dir, "remapping",
                   "Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv")
    # First run — produces baseline TSV
    subprocess.run(
        [
            "python", "scripts/benchmark_compare.py",
            "--manifest",   "backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv",
            "--remapped",   remapped_csv,
            "--assembly",   "equCab3",
            "--output-dir", str(run1),
        ],
        cwd="/home/tahmed/Equine80select_remapper",
        capture_output=True,
        check=True,
    )
    baseline_tsv = [str(run1 / f) for f in os.listdir(run1)
                    if f.endswith(".tsv") and "chrY" not in f and "chr0" not in f][0]

    run2 = tmp_path / "run2"
    run2.mkdir()
    # Second run — same data, with --baseline pointing at run1
    subprocess.run(
        [
            "python", "scripts/benchmark_compare.py",
            "--manifest",   "backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv",
            "--remapped",   remapped_csv,
            "--assembly",   "equCab3",
            "--output-dir", str(run2),
            "--baseline",   baseline_tsv,
        ],
        cwd="/home/tahmed/Equine80select_remapper",
        capture_output=True,
        check=True,
    )
    files = os.listdir(run2)
    assert any("_diff.txt" in f for f in files), f"_diff.txt not found in: {files}"


# ── load_remapped with new schema (anchor_ + tie_) ────────────────────────────

REMAPPED_CSV_NEW_SCHEMA = textwrap.dedent("""\
    Name,Chr_equCab3,MapInfo_equCab3,Strand_equCab3,anchor_equCab3,tie_equCab3
    SNP_auto,1,1000,+,topseq_n_probe,unique
    SNP_ambiguous,1,8000,+,topseq_n_probe,ambiguous
    SNP_topseq_only,2,5000,+,topseq_only,unique
    SNP_unmapped,0,0,N/A,N/A,N/A
    SNP_scaffold_resolved,3,7000,+,topseq_n_probe,scaffold_resolved
""")


@pytest.fixture
def remapped_file_new_schema(tmp_path):
    p = tmp_path / "remapped_new.csv"
    p.write_text(REMAPPED_CSV_NEW_SCHEMA)
    return str(p)


def test_load_remapped_new_schema_columns(remapped_file_new_schema):
    """New-schema CSV loads and produces the expected unified columns."""
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_new_schema, "equCab3")
    assert set(df.columns) >= {"Name", "remapped_chr", "remapped_pos",
                                "remapped_strand", "remapped_status"}
    # anchor_ and tie_ columns must be consumed, not left in the output
    assert "anchor_equCab3" not in df.columns
    assert "tie_equCab3" not in df.columns


def test_load_remapped_new_schema_status_mapped(remapped_file_new_schema):
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_new_schema, "equCab3")
    row = df.loc[df["Name"] == "SNP_auto"].iloc[0]
    assert row["remapped_status"] == "mapped"


def test_load_remapped_new_schema_status_ambiguous(remapped_file_new_schema):
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_new_schema, "equCab3")
    row = df.loc[df["Name"] == "SNP_ambiguous"].iloc[0]
    assert row["remapped_status"] == "ambiguous"


def test_load_remapped_new_schema_status_topseq_only(remapped_file_new_schema):
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_new_schema, "equCab3")
    row = df.loc[df["Name"] == "SNP_topseq_only"].iloc[0]
    assert row["remapped_status"] == "topseq_only"


def test_load_remapped_new_schema_status_unmapped(remapped_file_new_schema):
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_new_schema, "equCab3")
    row = df.loc[df["Name"] == "SNP_unmapped"].iloc[0]
    assert row["remapped_status"] == "unmapped"


def test_load_remapped_new_schema_scaffold_resolved_is_mapped(remapped_file_new_schema):
    """scaffold_resolved tie value with topseq_n_probe anchor → classified as 'mapped'."""
    from benchmark_compare import load_remapped
    df = load_remapped(remapped_file_new_schema, "equCab3")
    row = df.loc[df["Name"] == "SNP_scaffold_resolved"].iloc[0]
    assert row["remapped_status"] == "mapped"


def test_load_remapped_new_schema_ambiguous_classifies_correctly(
        remapped_file_new_schema, manifest_main_df):
    """End-to-end: new-schema ambiguous marker should classify as 'ambiguous'."""
    from benchmark_compare import load_remapped, compare_all
    # Add the ambiguous marker to a small manifest
    extra = pd.DataFrame([{
        "Name": "SNP_ambiguous", "manifest_chr": "1",
        "manifest_pos": 8000, "manifest_strand": "+",
    }])
    remapped_df = load_remapped(remapped_file_new_schema, "equCab3")
    result = compare_all(extra, remapped_df)
    assert result.loc[result["Name"] == "SNP_ambiguous", "result"].iloc[0] == "ambiguous"
