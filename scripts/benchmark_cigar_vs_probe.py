#!/usr/bin/env python3
"""
benchmark_cigar_vs_probe.py — Compare accuracy of three coordinate sources
against ground-truth manifest:

  1. CoordProbe_{assembly}        — probe-CIGAR coordinate (from probe alignment)
  2. Coord_TopSeqCIGAR_{assembly} — TopSeq-CIGAR coordinate (from TopGenomicSeq alignment)
  3. MapInfo_{assembly}           — final chosen coordinate
       (probe_cigar if CoordDelta < 2, topseq_cigar if CoordDelta ≥ 2)

Usage:
    python scripts/benchmark_cigar_vs_probe.py \\
        --manifest  backup_original/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \\
        --remapped  results_E80selv2_to_equCab3/..._remapped_equCab3.csv \\
        --assembly  equCab3
"""

import argparse
import contextlib
import io
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from benchmark_compare import (
    load_manifest,
    normalise_chr,
    classify_marker,
    CATEGORIES,
    _fmt,
    detect_assembly,
)
from remap_manifest import extract_candidates


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest",  required=True)
    p.add_argument("--remapped",  required=True)
    p.add_argument("--assembly",  default=None,
                   help="Assembly label (e.g. equCab3). Optional — auto-detected "
                        "from the remapped CSV header if omitted. Pass explicitly "
                        "only when the CSV contains multiple assemblies.")
    p.add_argument("--output-dir", default="./benchmark_out",
                   help="Output directory for the report file "
                        "(default: ./benchmark_out)")
    return p.parse_args()


def load_remapped_three_coords(path: str, assembly: str) -> pd.DataFrame:
    """
    Load the remapped CSV for a three-way probe/CIGAR/final coordinate comparison.

    Requires the current-schema columns `anchor_{assembly}` and `tie_{assembly}`.

    status derivation:
      anchor == "topseq_only"       → "topseq_only"
      tie    == "locus_unresolved"  → "locus_unresolved"
      anchor == "N/A"               → "unmapped"
      otherwise                     → "mapped"
    """
    col_chr         = f"Chr_{assembly}"
    col_pos_final   = f"MapInfo_{assembly}"
    col_pos_probe   = f"CoordProbe_{assembly}"
    col_pos_cigar   = f"Coord_TopSeqCIGAR_{assembly}"
    col_strand      = f"Strand_{assembly}"
    col_probe_strand = f"ProbeStrand_{assembly}"
    col_anchor      = f"anchor_{assembly}"
    col_tie         = f"tie_{assembly}"
    col_delta       = f"CoordDelta_{assembly}"
    col_source      = f"CoordSource_{assembly}"

    header = pd.read_csv(path, nrows=0)

    required_common = [col_chr, col_pos_final, col_pos_probe, col_pos_cigar,
                       col_strand, col_probe_strand, col_delta, col_source,
                       col_anchor, col_tie]
    for col in required_common:
        if col not in header.columns:
            raise ValueError(
                f"Missing expected column: {col!r}\n"
                f"Re-run remap_manifest.py to regenerate the remapped CSV with all "
                f"current-schema columns (including anchor_{assembly} + tie_{assembly})."
            )

    # TopGenomicSeq is needed to derive is_indel via extract_candidates;
    # match the pipeline's definition (len(AlleleA) != 1 or len(AlleleB) != 1
    # after "-" → "" normalisation).
    use_cols = ["Name", col_chr, col_pos_final, col_pos_probe, col_pos_cigar,
                col_strand, col_probe_strand, col_delta, col_source,
                col_anchor, col_tie]
    if "TopGenomicSeq" in header.columns:
        use_cols.append("TopGenomicSeq")

    df = pd.read_csv(
        path,
        dtype={col_chr: str},
        usecols=use_cols,
        low_memory=False,
    )

    # Build a unified status column
    def _derive_status(row):
        anchor = row[col_anchor]
        tie    = row[col_tie]
        # pandas reads "N/A" as NaN; treat NaN anchor as unmapped
        anchor_str = "" if pd.isna(anchor) else str(anchor)
        tie_str    = "" if pd.isna(tie)    else str(tie)
        if anchor_str == "topseq_only":
            return "topseq_only"
        if tie_str == "locus_unresolved":
            return "locus_unresolved"
        if anchor_str in ("N/A", ""):
            return "unmapped"
        return "mapped"
    df["status"] = df.apply(_derive_status, axis=1)
    df = df.drop(columns=[col_anchor, col_tie])

    # Derive is_indel using the same convention as remap_manifest.run_remapping
    # (extract_candidates normalises "-" → ""; then any allele with len != 1
    # is an indel).
    if "TopGenomicSeq" in df.columns:
        def _is_indel(ts):
            _, a, b, _ = extract_candidates(ts) if isinstance(ts, str) else (None, None, None, None)
            if a is None:
                return False
            return len(a) != 1 or len(b) != 1
        df["is_indel"] = df["TopGenomicSeq"].apply(_is_indel)
        df = df.drop(columns=["TopGenomicSeq"])
    else:
        df["is_indel"] = False

    df = df.rename(columns={
        col_chr:          "chr",
        col_pos_final:    "pos_final",
        col_pos_probe:    "pos_probe",
        col_pos_cigar:    "pos_cigar",
        col_strand:       "strand",
        col_probe_strand: "probe_strand",
        col_delta:        "coord_delta",
        col_source:       "coord_source",
    })
    df["chr"] = df["chr"].apply(normalise_chr)
    return df


def classify_with_pos(manifest_df, remapped_df, pos_col):
    """Classify markers using pos_col as the position.

    Threads `remapped_probe_strand` through to `classify_marker` so that the
    strand comparison uses probe alignment strand (which matches RefStrand ~97%)
    instead of TopSeq alignment strand (uncorrelated with RefStrand). Without
    this, `coord_correct_strand_wrong` would be reported on ~half of markers
    purely as a measurement artefact of comparing the wrong strand column.
    """
    merged = manifest_df.merge(remapped_df, on="Name", how="left")
    merged["chr"]          = merged["chr"].fillna("0")
    merged["strand"]       = merged["strand"].fillna("N/A")
    merged["probe_strand"] = merged["probe_strand"].fillna("N/A")
    merged["status"]       = merged["status"].fillna("unmapped")
    merged[pos_col]        = merged[pos_col].fillna(0)

    results = []
    for _, row in merged.iterrows():
        m = {
            "manifest_chr":    row["manifest_chr"],
            "manifest_pos":    row["manifest_pos"],
            "manifest_strand": row["manifest_strand"],
        }
        r = {
            "remapped_chr":          row["chr"],
            "remapped_pos":          row[pos_col],
            "remapped_strand":       row["strand"],
            "remapped_probe_strand": row["probe_strand"],
            "remapped_status":       row["status"],
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
    """Show accuracy of each coord source broken down by CoordDelta bucket.

    For each bucket prints three rows — total, SNV-only, indel-only — to make
    the per-source contribution to `final correct` reconcilable. The fifth
    column ('predicted') is derived per-marker from the `coord_source` chosen
    by the pipeline's selection cascade:
       coord_source == probe_cigar   → use probe result
       coord_source == topseq_cigar  → use CIGAR result
    Any drift between 'predicted' and 'final correct' attributes to
    post-decision adjustments (deletion minus-strand correction at
    remap_manifest.py:1729-1733, or final_pos shifts in determine_ref_alt_v2).
    """
    merged = main_df.merge(remapped_df, on="Name", how="left")
    merged["coord_delta"]  = merged["coord_delta"].fillna(-1).astype(int)
    merged["is_indel"]     = merged["is_indel"].fillna(False).astype(bool)
    merged["coord_source"] = merged["coord_source"].fillna("N/A")

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

    # "probe_cigar_indel_rescue" is a legacy label (pre-simplification) that
    # was unified into "probe_cigar". Accept it here so benchmarks can still
    # run against remapped CSVs generated before the rename. Safe to drop
    # once all in-use CSVs are regenerated.
    PROBE_SOURCES = {"probe_cigar", "probe_cigar_indel_rescue"}

    def _counts(names):
        n = len(names)
        if n == 0:
            return None
        sub = merged[merged["Name"].isin(names)]
        sources       = sub.set_index("Name")["coord_source"]
        probe_results = probe_by_name.reindex(names)
        cigar_results = cigar_by_name.reindex(names)
        final_results = final_by_name.reindex(names)
        # Predicted = probe-result for markers whose coord_source is probe-based,
        # cigar-result for the rest. Mirrors what the cascade rule selects.
        is_probe_src  = sources.reindex(names).isin(PROBE_SOURCES)
        predicted     = probe_results.where(is_probe_src, cigar_results)
        return {
            "n":         n,
            "probe":     (probe_results == "correct").sum(),
            "cigar":     (cigar_results == "correct").sum(),
            "final":     (final_results == "correct").sum(),
            "predicted": (predicted     == "correct").sum(),
        }

    col_w = 18
    label_w = 16
    print(f"  {'bucket':<{label_w}}  {'N':>6}"
          f"  {'probe correct':>{col_w}}"
          f"  {'CIGAR correct':>{col_w}}"
          f"  {'final correct':>{col_w}}"
          f"  {'predicted':>{col_w}}")
    print("  " + "-" * (label_w + 6 + (col_w + 2) * 4 + 6))

    def _row(label, c, indent=0):
        if c is None:
            return
        prefix = "  " * indent + label
        print(f"  {prefix:<{label_w}}  {c['n']:>6,}"
              f"  {_fmt(c['probe'],     c['n']):>{col_w}}"
              f"  {_fmt(c['cigar'],     c['n']):>{col_w}}"
              f"  {_fmt(c['final'],     c['n']):>{col_w}}"
              f"  {_fmt(c['predicted'], c['n']):>{col_w}}")

    for label, mask in buckets:
        all_names   = merged.loc[mask, "Name"]
        if len(all_names) == 0:
            continue
        snv_names   = merged.loc[mask & (~merged["is_indel"]), "Name"]
        indel_names = merged.loc[mask &   merged["is_indel"],  "Name"]

        _row(label,   _counts(all_names))
        _row("SNV",   _counts(snv_names),   indent=1)
        _row("indel", _counts(indel_names), indent=1)
    print()


def main():
    args = parse_args()

    print(f"[benchmark] Loading manifest: {args.manifest}")
    main_df, chry_df, chr0_df = load_manifest(args.manifest)
    print(f"[benchmark]   Benchmarked={len(main_df):,}  Chr=Y={len(chry_df):,}  Chr=0={len(chr0_df):,}")

    assembly = detect_assembly(args.remapped, args.assembly)
    if args.assembly is None:
        print(f"[benchmark] Auto-detected assembly: {assembly}")

    print(f"[benchmark] Loading remapped: {args.remapped}")
    remapped_df = load_remapped_three_coords(args.remapped, assembly)

    print("[benchmark] Classifying with probe coord (CoordProbe)...")
    probe_result = classify_with_pos(main_df, remapped_df, "pos_probe")

    print("[benchmark] Classifying with CIGAR coord (Coord_TopSeqCIGAR)...")
    cigar_result = classify_with_pos(main_df, remapped_df, "pos_cigar")

    print("[benchmark] Classifying with final coord (MapInfo)...")
    final_result = classify_with_pos(main_df, remapped_df, "pos_final")

    total = len(main_df)

    # Buffer the report sections so we can both write them to a file and
    # replay them to stdout (parity with benchmark_compare.py).
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
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
        print("  CoordDelta = |probe_cigar_coord − topseq_cigar_coord|; "
              "−1 = delta unavailable (topseq_only has no probe; probe_only has no TopSeq)")
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

        # Strategy comparison: accuracy under three coord-selection strategies
        print(f"{'='*96}")
        print("STRATEGY COMPARISON  (overall accuracy under each coord-selection strategy)")
        print(f"{'='*96}\n")
        strategies = [
            ("probe-only",       "always use CoordProbe",                          probe_result),
            ("topseq-only",      "always use Coord_TopSeqCIGAR",                   cigar_result),
            ("hybrid (current)", "probe_cigar if CoordDelta<2 else topseq_cigar",  final_result),
        ]
        for label, desc, df in strategies:
            n_correct = (df["result"] == "correct").sum()
            print(f"  {label:<18}  {desc:<48}  {_fmt(n_correct, total):>22}")
        print()

        probe_correct = (probe_result["result"] == "correct").sum()
        cigar_correct = (cigar_result["result"] == "correct").sum()
        final_correct = (final_result["result"] == "correct").sum()
        best_single   = max(probe_correct, cigar_correct)
        best_name     = "topseq-only" if cigar_correct >= probe_correct else "probe-only"
        gain          = final_correct - best_single
        gain_pp       = 100.0 * gain / total
        print(f"  Gain from hybrid over best single strategy ({best_name}): "
              f"{gain:+,} markers ({gain_pp:+.2f} pp)")
        print()

    # Write the report file
    os.makedirs(args.output_dir, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join(args.output_dir,
                               f"benchmark_cigar_vs_probe_{ts}_report.txt")
    with open(report_path, "w") as f:
        f.write(buf.getvalue())

    # Replay buffered content to the real stdout
    sys.stdout.write(buf.getvalue())
    print(f"[cigar_vs_probe] Report written to: {report_path}")


if __name__ == "__main__":
    main()
