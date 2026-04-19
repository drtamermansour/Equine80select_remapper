"""
prepare_benchmark_inputs.py — pre-stage inputs for the liftOver/CrossMap benchmark.

Given:
  - v1 Illumina manifest (mixed-assembly markers; the GenomeBuild=2 subset is our
    equCab2 benchmark set)
  - v2 Illumina manifest (all markers on equCab3; this is the ground truth)

Writes, under --output-dir:
  - v1_equCab2_subset.csv   a valid Illumina manifest containing only the
                            GenomeBuild=2 rows from v1. Suitable as input to
                            run_pipeline.sh.
  - v1_equCab2.bed          BED file of the equCab2 coordinates of those rows,
                            chr-prefixed for UCSC tools.
  - ground_truth.tsv        Name\\tchr_equCab3\\tpos_equCab3 for the markers that
                            also appear in v2 and have valid coordinates there.
"""

import argparse
import os
import sys
from typing import Iterable, List, Set, Tuple

import pandas as pd


# ── Manifest parsing ─────────────────────────────────────────────────────────

def split_manifest_sections(path: str) -> Tuple[List[str], str, List[str], List[str]]:
    """Split an Illumina manifest into four parts so we can filter-and-rewrite.

    Returns (preamble, column_header, data_lines, trailer) where
      preamble      = everything from the start through the ``[Assay]`` marker line,
      column_header = the ``IlmnID,Name,...`` line immediately after ``[Assay]``,
      data_lines    = the data rows between column_header and ``[Controls]``,
      trailer       = everything from ``[Controls]`` onwards.
    """
    with open(path) as f:
        lines = f.readlines()
    assay_idx = None
    controls_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip().rstrip(",")
        if stripped == "[Assay]" and assay_idx is None:
            assay_idx = i
        elif stripped == "[Controls]":
            controls_idx = i
            break
    if assay_idx is None:
        raise ValueError(f"No '[Assay]' section found in {path!r}")
    if controls_idx is None:
        controls_idx = len(lines)  # no [Controls] trailer; data runs to EOF
    preamble = lines[: assay_idx + 1]
    if assay_idx + 1 >= len(lines):
        raise ValueError(f"'[Assay]' marker is the last line of {path!r}")
    column_header = lines[assay_idx + 1]
    data_lines = lines[assay_idx + 2 : controls_idx]
    trailer = lines[controls_idx:]
    return preamble, column_header, data_lines, trailer


def read_manifest_df(path: str) -> pd.DataFrame:
    """Read an Illumina manifest's ``[Assay]`` data section into a DataFrame (all str)."""
    preamble, column_header, data_lines, _ = split_manifest_sections(path)
    # Write to an in-memory buffer and let pandas parse it uniformly.
    import io
    buf = io.StringIO(column_header + "".join(data_lines))
    return pd.read_csv(buf, dtype=str, keep_default_na=False)


# ── Filters ──────────────────────────────────────────────────────────────────

def _nonzero(series: pd.Series) -> pd.Series:
    """Boolean mask: value is a non-empty, non-zero string (interpreted as int)."""
    s = series.astype(str).str.strip()
    # empty / NaN / '0' / 'nan' fail
    return s.ne("") & s.ne("0") & s.str.lower().ne("nan")


def filter_genome_build_2(df: pd.DataFrame) -> pd.DataFrame:
    """Keep rows with GenomeBuild == '2' and valid Chr + MapInfo."""
    gb = df["GenomeBuild"].astype(str).str.strip()
    chr_ok = df["Chr"].notna() & _nonzero(df["Chr"])
    pos_ok = df["MapInfo"].notna() & _nonzero(df["MapInfo"])
    mask = gb.eq("2") & chr_ok & pos_ok
    return df.loc[mask].reset_index(drop=True)


def extract_ground_truth(v2_df: pd.DataFrame, names: Set[str]) -> pd.DataFrame:
    """Return v2 rows whose Name is in *names* and whose Chr+MapInfo are valid."""
    name_ok = v2_df["Name"].isin(names)
    chr_ok = v2_df["Chr"].notna() & _nonzero(v2_df["Chr"])
    pos_ok = v2_df["MapInfo"].notna() & _nonzero(v2_df["MapInfo"])
    return v2_df.loc[name_ok & chr_ok & pos_ok, ["Name", "Chr", "MapInfo"]].reset_index(drop=True)


# ── Writers ──────────────────────────────────────────────────────────────────

def df_to_bed_lines(df: pd.DataFrame) -> List[str]:
    """Build UCSC 0-based half-open BED records: ``chrN\\tpos-1\\tpos\\tName``."""
    out = []
    for _, row in df.iterrows():
        chrom = str(row["Chr"]).strip()
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        pos = int(row["MapInfo"])
        out.append(f"{chrom}\t{pos - 1}\t{pos}\t{row['Name']}\n")
    return out


def write_filtered_manifest(preamble: Iterable[str], column_header: str,
                            keep_names: Set[str], data_lines: Iterable[str],
                            trailer: Iterable[str], outpath: str) -> None:
    """Write a new Illumina manifest containing only the rows whose Name ∈ keep_names.

    The first column of the manifest is IlmnID, second is Name. We parse the
    column-header line to find the Name column index rather than assuming 2.
    """
    headers = column_header.rstrip("\r\n").split(",")
    try:
        name_idx = headers.index("Name")
    except ValueError:
        raise ValueError("Column 'Name' not found in manifest header")
    with open(outpath, "w") as f:
        f.writelines(preamble)
        f.write(column_header)
        for line in data_lines:
            parts = line.rstrip("\r\n").split(",")
            if len(parts) > name_idx and parts[name_idx] in keep_names:
                f.write(line)
        f.writelines(trailer)


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--v1", required=True, help="v1 (equCab2-designed) Illumina manifest CSV")
    p.add_argument("--v2", required=True, help="v2 (equCab3-native) Illumina manifest CSV (ground truth)")
    p.add_argument("-o", "--output-dir", required=True, help="Output directory (will be created)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # v1 — read, filter, preserve Illumina-manifest envelope for the subset.
    preamble_v1, header_v1, data_v1, trailer_v1 = split_manifest_sections(args.v1)
    v1_df = read_manifest_df(args.v1)
    print(f"[prep] v1 total data rows:          {len(v1_df):,}")

    v1_filtered = filter_genome_build_2(v1_df)
    print(f"[prep] v1 GenomeBuild=2 + valid Chr/MapInfo: {len(v1_filtered):,}")

    keep_names = set(v1_filtered["Name"].tolist())
    subset_path = os.path.join(args.output_dir, "v1_equCab2_subset.csv")
    write_filtered_manifest(preamble_v1, header_v1, keep_names, data_v1, trailer_v1, subset_path)
    print(f"[prep] wrote subset manifest:       {subset_path}")

    # BED file.
    bed_path = os.path.join(args.output_dir, "v1_equCab2.bed")
    with open(bed_path, "w") as f:
        f.writelines(df_to_bed_lines(v1_filtered))
    print(f"[prep] wrote equCab2 BED:           {bed_path}")

    # v2 — ground truth.
    v2_df = read_manifest_df(args.v2)
    print(f"[prep] v2 total data rows:          {len(v2_df):,}")
    gt = extract_ground_truth(v2_df, keep_names)
    gt_path = os.path.join(args.output_dir, "ground_truth.tsv")
    gt.to_csv(gt_path, sep="\t", index=False,
              header=["Name", "chr_equCab3", "pos_equCab3"])
    print(f"[prep] wrote ground truth:          {gt_path}")
    print(f"[prep] ground-truth rows:           {len(gt):,}  "
          f"(dropped {len(keep_names) - len(gt):,} markers missing from v2 or with Chr/MapInfo=0)")


if __name__ == "__main__":
    main()
