#!/usr/bin/env python3
"""
benchmark_cigar_vs_probe.py — Compare accuracy of three coordinate sources
against ground-truth manifest:

  1. CoordProbe_{assembly}        — raw probe-derived coordinate (before any override)
  2. Coord_TopSeqCIGAR_{assembly} — CIGAR-derived coordinate from TopGenomicSeq
  3. MapInfo_{assembly}           — final chosen coordinate (probe if delta<2, CIGAR otherwise)

Usage:
    python scripts/benchmark_cigar_vs_probe.py \\
        --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \\
        --remapped  results_E80selv2_to_equCab3/..._remapped_equCab3.csv \\
        --assembly  equCab3
"""

import argparse
import sys
import os

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from benchmark_compare import (
    load_manifest,
    normalise_chr,
    classify_marker,
    CATEGORIES,
    _fmt,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest",  required=True)
    p.add_argument("--remapped",  required=True)
    p.add_argument("--assembly",  required=True)
    return p.parse_args()


def load_remapped_three_coords(path: str, assembly: str) -> pd.DataFrame:
    col_chr         = f"Chr_{assembly}"
    col_pos_final   = f"MapInfo_{assembly}"
    col_pos_probe   = f"CoordProbe_{assembly}"
    col_pos_cigar   = f"Coord_TopSeqCIGAR_{assembly}"
    col_strand      = f"Strand_{assembly}"
    col_status      = f"MappingStatus_{assembly}"
    col_delta       = f"CoordDelta_{assembly}"
    col_source      = f"CoordSource_{assembly}"

    header = pd.read_csv(path, nrows=0)
    for col in [col_chr, col_pos_final, col_pos_probe, col_pos_cigar,
                col_strand, col_status, col_delta, col_source]:
        if col not in header.columns:
            raise ValueError(
                f"Missing expected column: {col!r}\n"
                f"Re-run remap_manifest.py to regenerate the remapped CSV with all columns."
            )

    df = pd.read_csv(
        path,
        dtype={col_chr: str},
        usecols=["Name", col_chr, col_pos_final, col_pos_probe, col_pos_cigar,
                 col_strand, col_status, col_delta, col_source],
        low_memory=False,
    )
    df = df.rename(columns={
        col_chr:       "chr",
        col_pos_final: "pos_final",
        col_pos_probe: "pos_probe",
        col_pos_cigar: "pos_cigar",
        col_strand:    "strand",
        col_status:    "status",
        col_delta:     "coord_delta",
        col_source:    "coord_source",
    })
    df["chr"] = df["chr"].apply(normalise_chr)
    return df


def classify_with_pos(manifest_df, remapped_df, pos_col):
    """Classify markers using pos_col as the position."""
    merged = manifest_df.merge(remapped_df, on="Name", how="left")
    merged["chr"]     = merged["chr"].fillna("0")
    merged["strand"]  = merged["strand"].fillna("N/A")
    merged["status"]  = merged["status"].fillna("unmapped")
    merged[pos_col]   = merged[pos_col].fillna(0)

    results = []
    for _, row in merged.iterrows():
        m = {
            "manifest_chr":    row["manifest_chr"],
            "manifest_pos":    row["manifest_pos"],
            "manifest_strand": row["manifest_strand"],
        }
        r = {
            "remapped_chr":    row["chr"],
            "remapped_pos":    row[pos_col],
            "remapped_strand": row["strand"],
            "remapped_status": row["status"],
        }
        # Treat pos=0 as unmapped for probe/CIGAR coords that are unavailable
        if row[pos_col] == 0:
            r["remapped_chr"]    = "0"
            r["remapped_strand"] = "N/A"
        results.append(classify_marker(m, r))
    merged["result"] = results
    return merged


def print_comparison(dfs_labels, total):
    labels = [label for _, label in dfs_labels]
    col_w  = 22

    header = f"  {'Category':<32}" + "".join(f"  {l:>{col_w}}" for l in labels)
    print(header)
    print("  " + "-" * (32 + (col_w + 2) * len(labels)))
    for cat in CATEGORIES:
        counts = [(_fmt((df["result"] == cat).sum(), total)) for df, _ in dfs_labels]
        print(f"  {cat:<32}" + "".join(f"  {c:>{col_w}}" for c in counts))
    print()


def print_delta_stratification(probe_df, cigar_df, final_df, remapped_df, main_df):
    """Show accuracy of each coord source broken down by CoordDelta bucket."""
    merged = main_df.merge(remapped_df, on="Name", how="left")
    merged["coord_delta"] = merged["coord_delta"].fillna(-1).astype(int)

    probe_by_name = probe_df.set_index("Name")["result"]
    cigar_by_name = cigar_df.set_index("Name")["result"]
    final_by_name = final_df.set_index("Name")["result"]

    delta = merged["coord_delta"]
    buckets = [
        ("delta = 0",    delta == 0),
        ("delta = 1",    delta == 1),
        ("delta = 2",    delta == 2),
        ("delta = 3",    delta == 3),
        ("delta = 4",    delta == 4),
        ("delta = 5",    delta == 5),
        ("delta = 6",    delta == 6),
        ("delta = 7",    delta == 7),
        ("delta = 8",    delta == 8),
        ("delta = 9",    delta == 9),
        ("delta = 10",   delta == 10),
        ("delta > 10",   delta > 10),
        ("delta = -1",   delta == -1),
    ]

    col_w = 22
    print(f"  {'bucket':<14}  {'N':>6}"
          f"  {'probe correct':>{col_w}}"
          f"  {'CIGAR correct':>{col_w}}"
          f"  {'final correct':>{col_w}}")
    print("  " + "-" * (14 + 6 + (col_w + 2) * 3 + 6))

    for label, mask in buckets:
        names = merged.loc[mask, "Name"]
        n = len(names)
        if n == 0:
            continue
        n_probe = (probe_by_name.reindex(names) == "correct").sum()
        n_cigar = (cigar_by_name.reindex(names) == "correct").sum()
        n_final = (final_by_name.reindex(names) == "correct").sum()
        print(f"  {label:<14}  {n:>6,}"
              f"  {_fmt(n_probe, n):>{col_w}}"
              f"  {_fmt(n_cigar, n):>{col_w}}"
              f"  {_fmt(n_final, n):>{col_w}}")
    print()


def main():
    args = parse_args()

    print(f"[benchmark] Loading manifest: {args.manifest}")
    main_df, chry_df, chr0_df = load_manifest(args.manifest)
    print(f"[benchmark]   Benchmarked={len(main_df):,}  Chr=Y={len(chry_df):,}  Chr=0={len(chr0_df):,}")

    print(f"[benchmark] Loading remapped: {args.remapped}")
    remapped_df = load_remapped_three_coords(args.remapped, args.assembly)

    print("[benchmark] Classifying with probe coord (CoordProbe)...")
    probe_result = classify_with_pos(main_df, remapped_df, "pos_probe")

    print("[benchmark] Classifying with CIGAR coord (Coord_TopSeqCIGAR)...")
    cigar_result = classify_with_pos(main_df, remapped_df, "pos_cigar")

    print("[benchmark] Classifying with final coord (MapInfo)...")
    final_result = classify_with_pos(main_df, remapped_df, "pos_final")

    total = len(main_df)

    print(f"\n{'='*96}")
    print(f"THREE-WAY COMPARISON  (N={total:,} benchmarked markers)")
    print(f"{'='*96}\n")
    print_comparison(
        [(probe_result, "probe (CoordProbe)"),
         (cigar_result, "CIGAR (TopSeqCIGAR)"),
         (final_result, "final (MapInfo)")],
        total,
    )

    print(f"{'='*96}")
    print("ACCURACY BY COORD_DELTA BUCKET  (correct count only)")
    print(f"  CoordDelta = |probe_coord - CIGAR_coord|; -1 = no probe (topseq_only or soft-clip)")
    print(f"{'='*96}\n")
    print_delta_stratification(probe_result, cigar_result, final_result, remapped_df, main_df)

    # Summary of CoordSource distribution
    print(f"{'='*96}")
    print("COORD SOURCE DISTRIBUTION  (how many markers used each source)")
    print(f"{'='*96}\n")
    merged = main_df.merge(remapped_df, on="Name", how="left")
    for src, count in merged["coord_source"].value_counts().items():
        print(f"  {src:<12}  {count:>8,}  ({100*count/total:.1f}%)")
    print()


if __name__ == "__main__":
    main()
