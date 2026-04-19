import os
import sys
import textwrap

import pandas as pd
import pytest

BENCHMARK_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "scripts", "benchmark")
)
if BENCHMARK_DIR not in sys.path:
    sys.path.insert(0, BENCHMARK_DIR)

from benchmark_vs_liftover import (
    classify_verdict,
    strip_chr_prefix,
    load_bed_predictions,
    load_allele_map_predictions,
    parse_ours_spec,
    summarize,
)


# ── classify_verdict ─────────────────────────────────────────────────────────

def test_correct_exact_match():
    assert classify_verdict(
        pred_chr="1", pred_pos=100, truth_chr="1", truth_pos=100, mapped=True
    ) == "correct"


def test_wrong_pos_le_10bp():
    assert classify_verdict("1", 105, "1", 100, True) == "wrong_pos_le_10bp"
    assert classify_verdict("1", 110, "1", 100, True) == "wrong_pos_le_10bp"


def test_wrong_pos_le_1kb():
    assert classify_verdict("1", 500, "1", 100, True) == "wrong_pos_le_1kb"
    assert classify_verdict("1", 1100, "1", 100, True) == "wrong_pos_le_1kb"


def test_wrong_pos_gt_1kb():
    assert classify_verdict("1", 50_000, "1", 100, True) == "wrong_pos_gt_1kb"


def test_wrong_chr():
    assert classify_verdict("2", 100, "1", 100, True) == "wrong_chr"


def test_unmapped():
    assert classify_verdict(None, None, "1", 100, False) == "unmapped"


def test_unmapped_ignores_pos_fields():
    # Even if pred_chr/pred_pos are present, mapped=False wins.
    assert classify_verdict("1", 100, "1", 100, False) == "unmapped"


# ── strip_chr_prefix ─────────────────────────────────────────────────────────

def test_strip_chr_prefix_removes_leading_chr():
    assert strip_chr_prefix("chr1") == "1"
    assert strip_chr_prefix("chrX") == "X"


def test_strip_chr_prefix_leaves_bare_chr_unchanged():
    assert strip_chr_prefix("1") == "1"
    assert strip_chr_prefix("MT") == "MT"


# ── load_bed_predictions ─────────────────────────────────────────────────────

def test_bed_predictions_parse_and_strip_prefix(tmp_path):
    bed = tmp_path / "lifted.bed"
    bed.write_text(
        "chr1\t99\t100\tsnp1\n"
        "chrX\t249\t250\tsnp2\n"
        "4\t50\t51\tsnp3\n"  # already bare
    )
    preds = load_bed_predictions(str(bed))
    assert preds["snp1"] == ("1", 100)
    assert preds["snp2"] == ("X", 250)
    assert preds["snp3"] == ("4", 51)


def test_bed_predictions_skip_comments_and_blanks(tmp_path):
    bed = tmp_path / "lifted.bed"
    bed.write_text(
        "# comment\n"
        "\n"
        "chr1\t99\t100\tsnp1\n"
    )
    preds = load_bed_predictions(str(bed))
    assert list(preds) == ["snp1"]


# ── load_allele_map_predictions ──────────────────────────────────────────────

ALLELE_MAP_TSV = textwrap.dedent("""\
    chr\tpos\tsnp_id\tmanifest_alleles\tgenomic_alleles\tmanifest_ref\tgenomic_ref\tdecision
    1\t1000\tsnpA\tA,G\tA,G\tA\tA\tas_is
    X\t5000\tsnpC\tC,T\tG,A\tC\tG\tcomplement
""")


def test_allele_map_loads_name_chr_pos(tmp_path):
    """Allele-map TSV → {Name: (chr, pos)}; only markers that survived QC appear."""
    p = tmp_path / "allele_map.tsv"
    p.write_text(ALLELE_MAP_TSV)
    preds = load_allele_map_predictions(str(p))
    assert preds == {"snpA": ("1", 1000), "snpC": ("X", 5000)}


def test_allele_map_missing_marker_treated_as_unmapped(tmp_path):
    """Markers absent from the allele-map are (by omission) unmapped for that preset."""
    p = tmp_path / "allele_map.tsv"
    p.write_text(ALLELE_MAP_TSV)
    preds = load_allele_map_predictions(str(p))
    assert "snpB" not in preds  # snpB was filtered out → caller treats as unmapped


# ── parse_ours_spec ──────────────────────────────────────────────────────────

def test_parse_ours_spec_splits_label_and_path():
    assert parse_ours_spec("default:/path/to/allele_map.tsv") == ("default", "/path/to/allele_map.tsv")


def test_parse_ours_spec_rejects_missing_colon():
    with pytest.raises(ValueError, match="must be of the form LABEL:PATH"):
        parse_ours_spec("just_a_path.tsv")


def test_parse_ours_spec_allows_colons_in_path():
    """Windows-style paths or URIs can contain colons; only the first is the separator."""
    assert parse_ours_spec("strict:C:/data/allele_map.tsv") == ("strict", "C:/data/allele_map.tsv")


# ── summarize ────────────────────────────────────────────────────────────────

def test_summarize_tallies_each_bucket():
    rows = [
        ("m1", "correct"),
        ("m2", "correct"),
        ("m3", "wrong_pos_le_10bp"),
        ("m4", "wrong_chr"),
        ("m5", "unmapped"),
    ]
    df = pd.DataFrame(rows, columns=["Name", "verdict"])
    counts = summarize(df, "verdict")
    assert counts["correct"] == 2
    assert counts["wrong_pos_le_10bp"] == 1
    assert counts["wrong_chr"] == 1
    assert counts["unmapped"] == 1
    # Buckets that didn't occur still appear as 0
    assert counts["wrong_pos_le_1kb"] == 0
    assert counts["wrong_pos_gt_1kb"] == 0
