import os
import sys

import pandas as pd
import pytest

# scripts/ is already on sys.path via conftest; benchmark/ is a subpackage location.
BENCHMARK_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "scripts", "benchmark")
)
if BENCHMARK_DIR not in sys.path:
    sys.path.insert(0, BENCHMARK_DIR)

from prepare_benchmark_inputs import (
    filter_genome_build_2,
    df_to_bed_lines,
    extract_ground_truth,
    split_manifest_sections,
)


# ── filter_genome_build_2 ────────────────────────────────────────────────────

def _v1(rows):
    return pd.DataFrame(
        rows, columns=["Name", "GenomeBuild", "Chr", "MapInfo"]
    )


def test_filter_keeps_genome_build_2():
    df = _v1([("a", "2", "1", "100"), ("b", "2", "X", "5000")])
    out = filter_genome_build_2(df)
    assert list(out["Name"]) == ["a", "b"]


def test_filter_drops_genome_build_3():
    df = _v1([("a", "2", "1", "100"), ("b", "3", "1", "200")])
    out = filter_genome_build_2(df)
    assert list(out["Name"]) == ["a"]


def test_filter_drops_chr_zero():
    df = _v1([("a", "2", "1", "100"), ("b", "2", "0", "200")])
    out = filter_genome_build_2(df)
    assert list(out["Name"]) == ["a"]


def test_filter_drops_chr_empty_or_nan():
    df = _v1([("a", "2", "1", "100"), ("b", "2", "", "200"), ("c", "2", None, "300")])
    out = filter_genome_build_2(df)
    assert list(out["Name"]) == ["a"]


def test_filter_drops_mapinfo_zero_or_empty():
    df = _v1([("a", "2", "1", "100"), ("b", "2", "1", "0"), ("c", "2", "1", "")])
    out = filter_genome_build_2(df)
    assert list(out["Name"]) == ["a"]


# ── df_to_bed_lines ──────────────────────────────────────────────────────────

def test_bed_lines_are_0_based_half_open_with_chr_prefix():
    df = pd.DataFrame([("snp1", "1", "100"), ("snp2", "X", "250")],
                      columns=["Name", "Chr", "MapInfo"])
    lines = df_to_bed_lines(df)
    assert lines == ["chr1\t99\t100\tsnp1\n", "chrX\t249\t250\tsnp2\n"]


def test_bed_lines_preserve_existing_chr_prefix():
    """If Chr already starts with 'chr', don't double-prefix."""
    df = pd.DataFrame([("snp1", "chr1", "100")],
                      columns=["Name", "Chr", "MapInfo"])
    lines = df_to_bed_lines(df)
    assert lines == ["chr1\t99\t100\tsnp1\n"]


# ── extract_ground_truth ─────────────────────────────────────────────────────

def _v2(rows):
    return pd.DataFrame(rows, columns=["Name", "Chr", "MapInfo"])


def test_ground_truth_only_keeps_names_in_input_set():
    df = _v2([("a", "1", "100"), ("b", "2", "200"), ("c", "3", "300")])
    out = extract_ground_truth(df, {"a", "c"})
    assert list(out["Name"]) == ["a", "c"]


def test_ground_truth_drops_rows_with_invalid_coords():
    df = _v2([("a", "1", "100"), ("b", "0", "200"), ("c", "1", "0"), ("d", "", "400")])
    out = extract_ground_truth(df, {"a", "b", "c", "d"})
    assert list(out["Name"]) == ["a"]


def test_ground_truth_returns_string_chr_and_int_pos():
    df = _v2([("a", "1", "100"), ("b", "X", "5000")])
    out = extract_ground_truth(df, {"a", "b"})
    assert out.dtypes["Chr"] == object
    # Pos must be integer-compatible for later comparison
    assert all(isinstance(int(v), int) for v in out["MapInfo"])


# ── split_manifest_sections ──────────────────────────────────────────────────

MINIMAL_MANIFEST = """[Header],,
Descriptor File Name,TestSet,,
Assay Format,Infinium,,
[Assay],,,
IlmnID,Name,GenomeBuild,Chr,MapInfo
id1,snp1,2,1,100
id2,snp2,2,X,500
id3,snp3,3,1,900
[Controls],,,
con1,red,,,
"""


def test_split_returns_four_parts(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text(MINIMAL_MANIFEST)
    preamble, header, data, trailer = split_manifest_sections(str(p))
    assert preamble[-1].rstrip("\r\n,") == "[Assay]"
    assert header.startswith("IlmnID,Name,")
    assert len(data) == 3
    assert trailer[0].rstrip("\r\n,") == "[Controls]"


def test_split_preamble_includes_header_metadata(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text(MINIMAL_MANIFEST)
    preamble, _, _, _ = split_manifest_sections(str(p))
    joined = "".join(preamble)
    assert "[Header]" in joined
    assert "Descriptor File Name" in joined
