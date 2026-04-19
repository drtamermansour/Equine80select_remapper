"""
benchmark_vs_liftover.py — multi-preset comparison of our remapper vs. liftOver vs. CrossMap.

Reads:
  --ground-truth                TSV (Name, chr_equCab3, pos_equCab3) produced by
                                prepare_benchmark_inputs.py
  --ours  LABEL:allele_map.tsv  our pipeline's per-preset allele-map. Repeatable
                                (once per preset you want to score).
  --liftover-lifted / --liftover-unmapped   BED files from `liftOver`
  --crossmap-lifted / --crossmap-unmapped   BED files from `CrossMap bed`

Writes under --output-dir:
  - three_way.tsv           per-marker verdicts for every preset + liftOver + CrossMap
  - benchmark_summary.txt   method × verdict table + sidebar + disagreement sample
"""

import argparse
import os
import sys
from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd


# ── Verdict taxonomy ─────────────────────────────────────────────────────────

VERDICT_ORDER = [
    "correct",
    "wrong_pos_le_10bp",
    "wrong_pos_le_1kb",
    "wrong_pos_gt_1kb",
    "wrong_chr",
    "unmapped",
]


def classify_verdict(pred_chr, pred_pos, truth_chr, truth_pos, mapped: bool) -> str:
    """Categorise a single marker's prediction against ground truth."""
    if not mapped:
        return "unmapped"
    if str(pred_chr) != str(truth_chr):
        return "wrong_chr"
    delta = abs(int(pred_pos) - int(truth_pos))
    if delta == 0:
        return "correct"
    if delta <= 10:
        return "wrong_pos_le_10bp"
    if delta <= 1000:
        return "wrong_pos_le_1kb"
    return "wrong_pos_gt_1kb"


# ── Prediction loaders ───────────────────────────────────────────────────────

def strip_chr_prefix(chrom: str) -> str:
    """Normalise UCSC-style chromosome names to bare form (``chr1`` → ``1``)."""
    s = str(chrom)
    return s[3:] if s.startswith("chr") else s


def load_bed_predictions(bed_path: str) -> Dict[str, Tuple[str, int]]:
    """Read a BED file (liftOver / CrossMap output) into {Name: (chr, pos_1based)}.

    BED is 0-based half-open; we report pos as ``end`` (the 1-based coordinate of
    the single-base marker). Comment lines and blanks are skipped.
    """
    preds: Dict[str, Tuple[str, int]] = {}
    if not os.path.exists(bed_path):
        return preds
    with open(bed_path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) < 4:
                continue
            chrom = strip_chr_prefix(parts[0])
            end = int(parts[2])
            name = parts[3]
            preds[name] = (chrom, end)
    return preds


def load_unmapped_names(bed_path: str) -> set:
    names = set()
    if not os.path.exists(bed_path):
        return names
    with open(bed_path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) >= 4:
                names.add(parts[3])
    return names


def load_allele_map_predictions(path: str) -> Dict[str, Tuple[str, int]]:
    """Read qc/{prefix}_allele_map_{assembly}.tsv into {Name: (chr, pos)}.

    Only markers that survived the QC cascade appear in this file — absence
    therefore means ``unmapped`` for that preset.
    """
    df = pd.read_csv(path, sep="\t", dtype={"chr": str})
    preds: Dict[str, Tuple[str, int]] = {}
    for _, row in df.iterrows():
        preds[str(row["snp_id"])] = (str(row["chr"]), int(row["pos"]))
    return preds


def parse_ours_spec(spec: str) -> Tuple[str, str]:
    """Split a ``LABEL:PATH`` argument into its two components.

    Only the **first** colon separates label from path — Windows-style paths
    like ``C:/foo/bar`` stay intact.
    """
    if ":" not in spec:
        raise ValueError(f"--ours value {spec!r} must be of the form LABEL:PATH")
    label, path = spec.split(":", 1)
    return label, path


# ── Summary ──────────────────────────────────────────────────────────────────

def summarize(df: pd.DataFrame, column: str) -> Counter:
    c = Counter(df[column])
    for v in VERDICT_ORDER:
        c.setdefault(v, 0)
    return c


def _format_table(method_counts, total_per_method):
    col_widths = {v: max(len(v), 9) for v in VERDICT_ORDER}
    label_w = max(len("method"), *(len(m) for m in method_counts))
    header = ["method".ljust(label_w), "total".rjust(8)] + [v.rjust(col_widths[v]) for v in VERDICT_ORDER]
    lines = ["  ".join(header), "  ".join("-" * len(h) for h in header)]
    for method, counts in method_counts.items():
        row = [method.ljust(label_w), str(total_per_method[method]).rjust(8)]
        for v in VERDICT_ORDER:
            row.append(str(counts[v]).rjust(col_widths[v]))
        lines.append("  ".join(row))
    return "\n".join(lines)


def _sidebar_from_allele_map(path: str) -> str:
    """Derive strand-flip / indel counts directly from the allele-map's ``decision`` column."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    n_total = len(df)
    decisions = df["decision"].fillna("").astype(str)
    n_complement = int(decisions.isin(["complement", "indel_complement"]).sum())
    n_indel      = int(decisions.str.startswith("indel_").sum())
    return (
        "Sidebar — features only our tool provides (from the 'default' preset's allele map)\n"
        f"  Total markers our tool placed:        {n_total:,}\n"
        f"  Strand-flip (decision=complement):    {n_complement:,}  (liftOver/CrossMap: coord only)\n"
        f"  Indel markers placed:                 {n_indel:,}       (liftOver/CrossMap: point coord only)\n"
    )


def _disagreement_sample(three_way: pd.DataFrame, ours_label: str, k: int = 20) -> str:
    """Top-k markers where ``ours_label``'s verdict disagrees with liftOver / CrossMap."""
    v_ours  = f"verdict_{ours_label}"
    c_ours  = f"{ours_label}_chr"
    p_ours  = f"{ours_label}_pos"
    dis = three_way[
        (three_way[v_ours] != three_way["verdict_liftover"]) |
        (three_way[v_ours] != three_way["verdict_crossmap"])
    ].copy()
    if dis.empty:
        return f"Qualitative disagreement sample ({ours_label} vs. liftOver/CrossMap): no disagreements.\n"
    ours_wins = dis[
        (dis[v_ours] == "correct")
        & ((dis["verdict_liftover"] != "correct") | (dis["verdict_crossmap"] != "correct"))
    ]
    ours_loses = dis[
        (dis[v_ours] != "correct")
        & ((dis["verdict_liftover"] == "correct") | (dis["verdict_crossmap"] == "correct"))
    ]
    sample = pd.concat([ours_wins.head(k), ours_loses.head(k)])
    cols = [
        "Name", "truth_chr", "truth_pos",
        c_ours, p_ours, v_ours,
        "lift_chr", "lift_pos", "verdict_liftover",
        "cross_chr", "cross_pos", "verdict_crossmap",
    ]
    return (
        f"Qualitative disagreement sample — top markers where '{ours_label}' differs from lift/cross\n"
        "-----------------------------------------------------------------------------------------\n"
        f"  ({ours_label}_wins={len(ours_wins):,};  {ours_label}_loses={len(ours_loses):,})\n\n"
        + sample[cols].to_string(index=False) + "\n"
    )


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ground-truth",     required=True)
    p.add_argument("--ours", action="append", required=True, metavar="LABEL:PATH",
                   help="Per-preset allele-map TSV. Repeatable. Example: "
                        "--ours default:ours/qc/v1_equCab2_subset_allele_map_equCab3.tsv "
                        "--ours strict:ours/qc_strict/v1_equCab2_subset_allele_map_equCab3.tsv")
    p.add_argument("--liftover-lifted",  required=True)
    p.add_argument("--liftover-unmapped", required=True)
    p.add_argument("--crossmap-lifted",  required=True)
    p.add_argument("--crossmap-unmapped", required=True)
    p.add_argument("-o", "--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    gt = pd.read_csv(args.ground_truth, sep="\t", dtype={"chr_equCab3": str})
    print(f"[bench] ground-truth markers:  {len(gt):,}")

    ours_specs: List[Tuple[str, str]] = [parse_ours_spec(s) for s in args.ours]
    ours_preds: Dict[str, Dict[str, Tuple[str, int]]] = {}
    for label, path in ours_specs:
        ours_preds[label] = load_allele_map_predictions(path)
        print(f"[bench] ours[{label}]  placed markers: {len(ours_preds[label]):,}")

    lift_mapped   = load_bed_predictions(args.liftover_lifted)
    lift_unmapped = load_unmapped_names(args.liftover_unmapped)
    cross_mapped  = load_bed_predictions(args.crossmap_lifted)
    cross_unmapped = load_unmapped_names(args.crossmap_unmapped)
    print(f"[bench] liftOver  mapped: {len(lift_mapped):,};  unmapped: {len(lift_unmapped):,}")
    print(f"[bench] CrossMap  mapped: {len(cross_mapped):,};  unmapped: {len(cross_unmapped):,}")

    # Build the wide per-marker verdict frame.
    rows = []
    base_cols = ["Name", "truth_chr", "truth_pos"]
    ours_cols: List[str] = []
    for label, _ in ours_specs:
        ours_cols += [f"{label}_chr", f"{label}_pos", f"verdict_{label}"]
    tail_cols = [
        "lift_chr", "lift_pos", "verdict_liftover",
        "cross_chr", "cross_pos", "verdict_crossmap",
    ]

    for _, row in gt.iterrows():
        name = row["Name"]
        truth_chr, truth_pos = row["chr_equCab3"], int(row["pos_equCab3"])
        record = [name, truth_chr, truth_pos]

        # Each ours-preset
        for label, _ in ours_specs:
            pred = ours_preds[label].get(name)
            if pred is None:
                record += ["", "", "unmapped"]
            else:
                c, p = pred
                v = classify_verdict(c, p, truth_chr, truth_pos, True)
                record += [c, p, v]

        # liftOver
        if name in lift_mapped:
            lc, lp = lift_mapped[name]
            record += [lc, lp, classify_verdict(lc, lp, truth_chr, truth_pos, True)]
        else:
            record += ["", "", "unmapped"]

        # CrossMap
        if name in cross_mapped:
            cc, cp = cross_mapped[name]
            record += [cc, cp, classify_verdict(cc, cp, truth_chr, truth_pos, True)]
        else:
            record += ["", "", "unmapped"]

        rows.append(record)

    three_way = pd.DataFrame(rows, columns=base_cols + ours_cols + tail_cols)
    out_tsv = os.path.join(args.output_dir, "three_way.tsv")
    three_way.to_csv(out_tsv, sep="\t", index=False)
    print(f"[bench] wrote {out_tsv}")

    # Summary: one row per ours-preset, then liftOver, then CrossMap.
    method_counts = {}
    for label, _ in ours_specs:
        method_counts[f"ours ({label})"] = summarize(three_way, f"verdict_{label}")
    method_counts["liftOver"] = summarize(three_way, "verdict_liftover")
    method_counts["CrossMap"] = summarize(three_way, "verdict_crossmap")
    totals = {m: len(three_way) for m in method_counts}

    # Sidebar: use the "default" preset if present, otherwise the first.
    sidebar_label, sidebar_path = next(
        ((l, p) for l, p in ours_specs if l == "default"),
        ours_specs[0],
    )
    sidebar = _sidebar_from_allele_map(sidebar_path)

    summary_path = os.path.join(args.output_dir, "benchmark_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Benchmark: our tool vs. liftOver vs. CrossMap (equCab2 → equCab3)\n")
        f.write(f"Ground-truth markers: {len(three_way):,}\n")
        f.write(f"'ours' presets:       {', '.join(l for l, _ in ours_specs)}\n")
        f.write("=" * 88 + "\n\n")
        f.write("Head-to-head (coordinates only)\n")
        f.write("-" * 88 + "\n")
        f.write(_format_table(method_counts, totals) + "\n\n")
        f.write(sidebar + "\n")
        f.write("=" * 88 + "\n")
        f.write(_disagreement_sample(three_way, sidebar_label) + "\n")
    print(f"[bench] wrote {summary_path}")


if __name__ == "__main__":
    main()
