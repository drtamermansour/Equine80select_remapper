import os
import pytest
import pandas as pd
from benchmark_compare import (
    normalise_chr, classify_marker, write_diff,
    stratify_by_coord_delta,
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
            "--remapped",   os.path.join(results_dir,
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
            "--remapped",   os.path.join(results_dir,
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
    remapped_csv = os.path.join(results_dir,
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
