#!/usr/bin/env python3
"""
benchmark_compare.py — Compare remap_manifest.py output against a ground-truth
EquCab3-native Illumina manifest. Produces per-marker TSVs and a summary report.

Usage:
    python scripts/benchmark_compare.py \
        --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
        --remapped  results_E80selv2_to_equCab3/..._remapped_equCab3.csv \
        --assembly  equCab3 \
        [--output-dir results_E80selv2_to_equCab3/benchmark/] \
        [--baseline  results_E80selv2_to_equCab3/benchmark/benchmark_TIMESTAMP.tsv]
"""

import argparse
import os
import sys
from datetime import datetime

import pandas as pd

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Benchmark remap_manifest.py output against ground-truth manifest.")
    p.add_argument("--manifest",    required=True, help="Illumina manifest CSV (EquCab3-native, GenomeBuild=3)")
    p.add_argument("--remapped",    required=True, help="*_remapped_{assembly}.csv from remap_manifest.py")
    p.add_argument("--assembly",    required=True, help="Assembly label (e.g. equCab3)")
    p.add_argument("--output-dir",  default="./benchmark_out", help="Output directory (default: ./benchmark_out)")
    p.add_argument("--baseline",    default=None,  help="Path to prior benchmark_<timestamp>.tsv for diff")
    return p.parse_args()


# ── CHROMOSOME NORMALISATION ──────────────────────────────────────────────────

def normalise_chr(chrom: str) -> str:
    """Normalise X chromosome aliases to 'X'. All other values pass through unchanged."""
    if str(chrom).startswith("X_"):
        return "X"
    return str(chrom)


# ── COMPARISON ────────────────────────────────────────────────────────────────

def classify_marker(manifest_row: dict, remapped_row: dict) -> str:
    """
    Classify a single marker into one of 6 outcome categories.

    Category priority (highest first):
      unmapped  — remapped Chr=0 or Strand=N/A
      ambiguous — MappingStatus=ambiguous
      wrong_chr — Chr mismatch
      coord_off — Chr match but MapInfo differs
      coord_correct_strand_wrong — Chr+MapInfo match but Strand differs
      correct   — all three match
    """
    m_chr    = normalise_chr(manifest_row["manifest_chr"])
    m_pos    = manifest_row["manifest_pos"]
    m_strand = manifest_row["manifest_strand"]

    r_chr    = normalise_chr(remapped_row["remapped_chr"])
    r_pos    = remapped_row["remapped_pos"]
    r_strand = remapped_row["remapped_strand"]
    r_status = remapped_row["remapped_status"]

    # Priority 1: unmapped
    if r_chr == "0" or r_strand == "N/A":
        return "unmapped"

    # Priority 2: ambiguous (note: ambiguous short-circuits all coordinate checks,
    # so an ambiguous marker with a wrong chromosome is still counted as 'ambiguous',
    # not 'wrong_chr')
    if r_status == "ambiguous":
        return "ambiguous"

    # Priority 3: wrong chromosome
    if r_chr != m_chr:
        return "wrong_chr"

    # Priority 4: coordinate off
    try:
        if int(r_pos) != int(m_pos):
            return "coord_off"
    except (ValueError, TypeError):
        return "coord_off"

    # Priority 5: strand wrong
    if str(r_strand) != str(m_strand):
        return "coord_correct_strand_wrong"

    return "correct"


# ── MANIFEST LOADING ──────────────────────────────────────────────────────────

def _locate_assay_section(path: str) -> tuple[int, int | None]:
    """
    Returns (header_line, footer_line) as 0-based line indices.
    header_line  : the line containing the column names (line after [Assay])
    footer_line  : the line containing [Controls], or None if absent
    """
    header_line = None
    footer_line = None
    with open(path) as f:
        for i, line in enumerate(f):
            stripped = line.strip().rstrip(",")
            if stripped == "[Assay]":
                header_line = i + 1
            if stripped == "[Controls]":
                footer_line = i
                break
    if header_line is None:
        raise ValueError(f"[Assay] section not found in {path}")
    return header_line, footer_line


def load_manifest(path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load an Illumina manifest CSV and partition markers into three DataFrames:
      main_df  — autosome or X (benchmarked)
      chry_df  — Chr == Y
      chr0_df  — Chr == 0 or MapInfo == 0

    All DataFrames have columns: Name, manifest_chr, manifest_pos, manifest_strand
    Chr is normalised (X_NC_009175.3 → X) in all three.
    """
    header_line, footer_line = _locate_assay_section(path)
    nrows = (footer_line - header_line - 1) if footer_line is not None else None

    df = pd.read_csv(
        path,
        skiprows=header_line,
        nrows=nrows,
        dtype={"Chr": str, "MapInfo": str},
        usecols=["Name", "Chr", "MapInfo", "SourceStrand"],
        low_memory=False,
    )
    df = df.rename(columns={"Chr": "manifest_chr", "MapInfo": "manifest_pos",
                             "SourceStrand": "manifest_strand"})
    df["manifest_chr"] = df["manifest_chr"].apply(normalise_chr)
    df["manifest_pos"] = pd.to_numeric(df["manifest_pos"], errors="coerce").fillna(0).astype(int)
    # Normalise Illumina strand notation to +/-
    _ilmn_map = {"TOP": "+", "PLUS": "+", "BOT": "-", "MINUS": "-"}
    df["manifest_strand"] = df["manifest_strand"].str.upper().map(
        lambda s: _ilmn_map.get(s, s)
    )

    # Partition
    is_y    = df["manifest_chr"] == "Y"
    is_chr0 = (df["manifest_chr"] == "0") | (df["manifest_pos"] == 0)
    is_main = ~is_y & ~is_chr0

    return df[is_main].reset_index(drop=True), \
           df[is_y].reset_index(drop=True), \
           df[is_chr0].reset_index(drop=True)


# ── REMAPPED LOADING ──────────────────────────────────────────────────────────

def load_remapped(path: str, assembly: str) -> pd.DataFrame:
    """
    Load the remapped CSV produced by remap_manifest.py.
    Returns a DataFrame with columns:
      Name, remapped_chr, remapped_pos, remapped_strand, remapped_status
    and optionally coord_delta (if CoordDelta_{assembly} is present in the CSV).
    Chr is normalised (X aliases → X).

    Supports both the new column schema (anchor_{assembly} + tie_{assembly}) and the
    old schema (MappingStatus_{assembly}) for backward compatibility.  When both are
    present the new schema takes precedence.

    New schema → remapped_status mapping:
      anchor == "topseq_only"              → "topseq_only"
      tie    == "ambiguous"                → "ambiguous"
      anchor == "N/A" (unmapped)           → "unmapped"
      otherwise                            → "mapped"
    """
    col_chr    = f"Chr_{assembly}"
    col_pos    = f"MapInfo_{assembly}"
    col_strand = f"Strand_{assembly}"
    col_status = f"MappingStatus_{assembly}"   # old schema
    col_anchor = f"anchor_{assembly}"          # new schema
    col_tie    = f"tie_{assembly}"             # new schema
    col_delta  = f"CoordDelta_{assembly}"

    header = pd.read_csv(path, nrows=0)

    has_new_schema = col_anchor in header.columns and col_tie in header.columns
    has_old_schema = col_status in header.columns

    # Validate required columns
    required = ["Name", col_chr, col_pos, col_strand]
    if not has_new_schema and not has_old_schema:
        required += [col_status]  # will produce a meaningful missing-column error
    missing = [c for c in required if c not in header.columns]
    if missing:
        raise ValueError(
            f"Remapped CSV {path!r} is missing expected columns: {missing}\n"
            f"Check that --assembly matches the assembly used when running remap_manifest.py."
        )
    has_delta = col_delta in header.columns

    usecols = ["Name", col_chr, col_pos, col_strand]
    if has_new_schema:
        usecols += [col_anchor, col_tie]
    elif has_old_schema:
        usecols.append(col_status)
    if has_delta:
        usecols.append(col_delta)

    df = pd.read_csv(
        path,
        dtype={col_chr: str},
        usecols=usecols,
        low_memory=False,
    )

    # Build a unified remapped_status column
    if has_new_schema:
        def _derive_status(row):
            anchor = row[col_anchor]
            tie    = row[col_tie]
            # pandas reads "N/A" as NaN; treat NaN anchor as unmapped
            anchor_str = "" if pd.isna(anchor) else str(anchor)
            tie_str    = "" if pd.isna(tie)    else str(tie)
            if anchor_str == "topseq_only":
                return "topseq_only"
            if tie_str == "ambiguous":
                return "ambiguous"
            if anchor_str in ("N/A", ""):
                return "unmapped"
            return "mapped"
        df["remapped_status"] = df.apply(_derive_status, axis=1)
        df = df.drop(columns=[col_anchor, col_tie])
    else:
        df = df.rename(columns={col_status: "remapped_status"})

    rename_map = {
        col_chr:    "remapped_chr",
        col_pos:    "remapped_pos",
        col_strand: "remapped_strand",
    }
    if has_delta:
        rename_map[col_delta] = "coord_delta"
    df = df.rename(columns=rename_map)
    df["remapped_chr"] = df["remapped_chr"].apply(normalise_chr)
    df["remapped_pos"] = pd.to_numeric(df["remapped_pos"], errors="coerce")
    return df


# ── COMPARISON ────────────────────────────────────────────────────────────────

def compare_all(manifest_df: pd.DataFrame, remapped_df: pd.DataFrame) -> pd.DataFrame:
    """
    Left-join manifest onto remapped, classify each marker, compute coord_offset.
    Markers present in manifest but absent from remapped are classified as 'unmapped'.

    Returns a DataFrame with all manifest columns plus:
      remapped_chr, remapped_pos, remapped_strand, remapped_status, result, coord_offset
    """
    merged = manifest_df.merge(remapped_df, on="Name", how="left")

    # Fill missing remapped rows (not in remapped CSV at all)
    merged["remapped_chr"]    = merged["remapped_chr"].fillna("0")
    merged["remapped_pos"]    = merged["remapped_pos"].fillna(0)
    merged["remapped_strand"] = merged["remapped_strand"].fillna("N/A")
    merged["remapped_status"] = merged["remapped_status"].fillna("unmapped")

    results = []
    offsets = []
    for _, row in merged.iterrows():
        m = {k: row[k] for k in ("manifest_chr", "manifest_pos", "manifest_strand")}
        r = {k: row[k] for k in ("remapped_chr", "remapped_pos",
                                  "remapped_strand", "remapped_status")}
        cat = classify_marker(m, r)
        results.append(cat)
        if cat == "coord_off":
            try:
                offsets.append(int(row["remapped_pos"]) - int(row["manifest_pos"]))
            except (ValueError, TypeError):
                offsets.append(None)
        else:
            offsets.append(None)

    merged["result"]       = results
    merged["coord_offset"] = offsets
    if "coord_delta" in merged.columns:
        merged["coord_delta"] = merged["coord_delta"].fillna(-1).astype(int)
    return merged


# ── OUTPUT ────────────────────────────────────────────────────────────────────

_TSV_COLS = [
    "Name",
    "manifest_chr", "manifest_pos", "manifest_strand",
    "remapped_chr",  "remapped_pos",  "remapped_strand", "remapped_status",
    "result", "coord_offset", "coord_delta",
]

CATEGORIES = [
    "correct",
    "coord_correct_strand_wrong",
    "coord_off",
    "wrong_chr",
    "unmapped",
    "ambiguous",
]


def _compute_transitions(
    curr: pd.Series, base: pd.Series
) -> tuple[int, dict[tuple[str, str], int]]:
    """Return (changed_count, transitions_dict) between two result series indexed by Name."""
    common = curr.index.intersection(base.index)
    changed = common[curr[common] != base[common]]
    transitions: dict[tuple[str, str], int] = {}
    for name in changed:
        key = (base[name], curr[name])
        transitions[key] = transitions.get(key, 0) + 1
    return len(changed), transitions


def write_tsv(df: pd.DataFrame, path: str) -> None:
    """Write result DataFrame to a tab-separated file using the standard column order."""
    cols = [c for c in _TSV_COLS if c in df.columns]
    df[cols].to_csv(path, sep="\t", index=False)


def load_baseline(path: str) -> pd.DataFrame:
    """Load a previously written benchmark TSV for diff comparison."""
    return pd.read_csv(path, sep="\t", dtype={"manifest_chr": str, "remapped_chr": str})


def _fmt(n: int, total: int) -> str:
    pct = 100 * n / total if total else 0.0
    return f"{n:>10,}  ({pct:5.1f}%)"


def _marker_type(name: str) -> str:
    if str(name).startswith("Affx-"):
        return "AFFX"
    if "ilmndup" in str(name):
        return "ilmndup"
    return "standard"


def stratify_by_coord_delta(result_df: pd.DataFrame) -> list[dict] | None:
    """
    Stratify accuracy by CoordDelta bucket.

    CoordDelta = |probe_coord - CIGAR_coord|; -1 means CIGAR coord was unavailable
    (SNP target fell in a soft-clipped region).

    Returns a list of dicts (one per non-empty bucket), each with:
      label, n, coord_accurate (correct + coord_correct_strand_wrong), correct

    Returns None if the coord_delta column is absent (old remapped CSVs).
    """
    if "coord_delta" not in result_df.columns:
        return None

    delta = result_df["coord_delta"]
    buckets = [
        ("delta = 0",    delta == 0),
        ("delta = 1",    delta == 1),
        ("delta = 2-10", (delta >= 2) & (delta <= 10)),
        ("delta > 10",   delta > 10),
        ("delta = -1",   delta == -1),
    ]

    rows = []
    for label, mask in buckets:
        subset = result_df[mask]
        n = len(subset)
        if n == 0:
            continue
        n_correct       = int((subset["result"] == "correct").sum())
        n_strand_wrong  = int((subset["result"] == "coord_correct_strand_wrong").sum())
        rows.append({
            "label":          label,
            "n":              n,
            "coord_accurate": n_correct + n_strand_wrong,
            "correct":        n_correct,
        })
    return rows


def write_report(
    result_df: pd.DataFrame,
    chry_df: pd.DataFrame,
    chr0_df: pd.DataFrame,
    path: str,
    assembly: str,
    manifest_path: str,
    remapped_path: str,
    baseline_df: pd.DataFrame | None = None,
) -> None:
    """Write the human-readable benchmark report."""
    total_manifest = len(result_df) + len(chry_df) + len(chr0_df)
    benchmarked    = len(result_df)
    lines = []

    def w(s=""):
        lines.append(s)

    w(f"Benchmark Report — assembly: {assembly}")
    w(f"Run: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    w(f"Manifest:  {manifest_path}")
    w(f"Remapped:  {remapped_path}")
    w("-" * 60)
    w()
    w("SCOPE")
    w(f"  Total manifest markers:        {total_manifest:>10,}")
    w(f"  Benchmarked (autosome + X):    {benchmarked:>10,}")
    w(f"  Excluded Chr=Y:                {len(chry_df):>10,}")
    w(f"  Excluded Chr=0:                {len(chr0_df):>10,}")
    w()
    w("HEADLINE COUNTS")
    for cat in CATEGORIES:
        n = (result_df["result"] == cat).sum()
        w(f"  {cat:<32} {_fmt(n, benchmarked)}")
    w()
    w("BREAKDOWN BY MARKER TYPE")
    for mtype in ("standard", "AFFX", "ilmndup"):
        subset = result_df[result_df["Name"].apply(_marker_type) == mtype]
        if len(subset) == 0:
            continue
        w(f"  {mtype} ({len(subset):,} markers):")
        for cat in CATEGORIES:
            n = (subset["result"] == cat).sum()
            if n > 0:
                w(f"    {cat:<30} {_fmt(n, len(subset))}")
    w()

    # Coordinate offset distribution (coord_off markers only)
    coord_off = result_df[result_df["result"] == "coord_off"].copy()
    if len(coord_off) > 0:
        w("COORDINATE OFFSET DISTRIBUTION  (coord_off markers only)")
        offsets = coord_off["coord_offset"].dropna().abs()
        buckets = [
            ("offset = 1 bp",                      offsets == 1),
            ("offset = 2–10 bp",                   (offsets >= 2) & (offsets <= 10)),
            ("offset = 11–50 bp",                  (offsets >= 11) & (offsets <= 50)),
            ("offset = 51 bp (probe-strand bug)",  offsets == 51),
            ("offset = 52+ bp",                    offsets >= 52),
        ]
        for label, mask in buckets:
            w(f"  {label:<40} {mask.sum():>8,}")
        w()

    # CoordDelta stratification (only when column present)
    delta_rows = stratify_by_coord_delta(result_df)
    if delta_rows is not None:
        w("ACCURACY STRATIFIED BY COORD_DELTA")
        w("  CoordDelta = |probe_coord - CIGAR_coord|; -1 = CIGAR target in soft clip")
        w("  coord-accurate = correct + coord_correct_strand_wrong")
        w()
        w(f"  {'bucket':<14}  {'N':>8}  {'coord-accurate':>20}  {'correct':>20}")
        w("  " + "-" * 68)
        for r in delta_rows:
            ca_str = _fmt(r["coord_accurate"], r["n"])
            co_str = _fmt(r["correct"],        r["n"])
            w(f"  {r['label']:<14}  {r['n']:>8,}  {ca_str}  {co_str}")
        w()

    # Diff section
    if baseline_df is not None:
        w("-" * 60)
        w("DIFF VS BASELINE")
        curr = result_df.set_index("Name")["result"]
        base = baseline_df.set_index("Name")["result"]
        n_changed, transitions = _compute_transitions(curr, base)
        w(f"  Markers that changed category: {n_changed:>8,}")
        w()
        for (from_cat, to_cat), count in sorted(transitions.items(), key=lambda x: -x[1]):
            w(f"  {from_cat:<32} → {to_cat:<32} {count:>8,}")
        w()

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_diff(result_df: pd.DataFrame, baseline_df: pd.DataFrame, path: str) -> None:
    """Write a standalone diff file listing category transitions vs baseline."""
    curr = result_df.set_index("Name")["result"]
    base = baseline_df.set_index("Name")["result"]
    n_changed, transitions = _compute_transitions(curr, base)

    lines = []
    lines.append(f"Markers that changed category: {n_changed:,}")
    lines.append("")
    for (from_cat, to_cat), count in sorted(transitions.items(), key=lambda x: -x[1]):
        lines.append(f"  {from_cat:<32} → {to_cat:<32} {count:>8,}")
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    print(f"[benchmark] Loading manifest: {args.manifest}")
    main_df, chry_df, chr0_df = load_manifest(args.manifest)
    print(f"[benchmark]   Benchmarked={len(main_df):,}  Chr=Y={len(chry_df):,}  Chr=0={len(chr0_df):,}")

    print(f"[benchmark] Loading remapped: {args.remapped}")
    remapped_df = load_remapped(args.remapped, args.assembly)

    print("[benchmark] Classifying markers...")
    result_df = compare_all(main_df, remapped_df)

    # Classify chrY and chr0 rows for their side TSVs
    chry_result = compare_all(chry_df, remapped_df)
    chry_result["result"] = "chrY"
    chr0_result = compare_all(chr0_df, remapped_df)
    chr0_result["result"] = "chr0"

    # Write TSVs
    main_tsv = os.path.join(args.output_dir, f"benchmark_{ts}.tsv")
    chry_tsv = os.path.join(args.output_dir, f"benchmark_{ts}_chrY.tsv")
    chr0_tsv = os.path.join(args.output_dir, f"benchmark_{ts}_chr0.tsv")
    write_tsv(result_df,   main_tsv)
    write_tsv(chry_result, chry_tsv)
    write_tsv(chr0_result, chr0_tsv)
    print(f"[benchmark] TSVs written to: {args.output_dir}")

    # Load baseline if provided
    baseline_df = None
    if args.baseline:
        print(f"[benchmark] Loading baseline: {args.baseline}")
        baseline_df = load_baseline(args.baseline)

    # Write report
    report_path = os.path.join(args.output_dir, f"benchmark_{ts}_report.txt")
    write_report(
        result_df, chry_result, chr0_result,
        report_path,
        assembly=args.assembly,
        manifest_path=args.manifest,
        remapped_path=args.remapped,
        baseline_df=baseline_df,
    )
    print(f"[benchmark] Report: {report_path}")

    # Write standalone diff file if baseline was provided
    if baseline_df is not None:
        diff_path = os.path.join(args.output_dir, f"benchmark_{ts}_diff.txt")
        write_diff(result_df, baseline_df, diff_path)
        print(f"[benchmark] Diff:   {diff_path}")

    # Print report to stdout
    print()
    with open(report_path) as f:
        print(f.read())


if __name__ == "__main__":
    main()
