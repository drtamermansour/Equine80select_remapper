"""
filter_scaffold_haplotypes.py — Filter scaffold_summary.tsv to candidate alt haplotypes.

Applies alignment-based thresholds to the output of scaffold_haplotype_analyzer.py and
writes a filtered TSV of scaffolds that are likely alternative haplotypes of placed chromosomes.

Default thresholds are Tier 1 (high confidence). All thresholds are adjustable via CLI.
See docs/scaffold_haplotype_thresholds.md for the full classification rationale.

Usage:
  python scripts/filter_scaffold_haplotypes.py \\
      -i remap_assessment/scaffold_haplotype_analysis/scaffold_summary.tsv \\
      -o remap_assessment/scaffold_haplotype_analysis/alt_haplotype_candidates.tsv
"""

import argparse
import os
import sys

import pandas as pd


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Filter scaffold_summary.tsv to likely alt-haplotype candidates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input",  required=True,
                   help="scaffold_summary.tsv from scaffold_haplotype_analyzer.py")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV of passing scaffolds")

    # Tier 1 thresholds as defaults
    p.add_argument("--min-identity",      type=float, default=99.0,
                   help="Minimum identity_pct")
    p.add_argument("--min-query-cov",     type=float, default=80.0,
                   help="Minimum query_coverage_pct")
    p.add_argument("--max-span-ratio",    type=float, default=1.5,
                   help="Maximum span_to_scaffold_ratio")
    p.add_argument("--min-mapq",          type=int,   default=40,
                   help="Minimum max_mapq")
    p.add_argument("--max-blocks",        type=int,   default=5,
                   help="Maximum n_alignment_blocks")
    return p.parse_args()


# ── FILTER ────────────────────────────────────────────────────────────────────

def filter_scaffolds(df, args):
    """Apply threshold filters; return filtered DataFrame and per-filter counts."""
    counts = {"input": len(df)}

    df = df[df["identity_pct"]          >= args.min_identity];   counts["identity"]   = len(df)
    df = df[df["query_coverage_pct"]    >= args.min_query_cov];  counts["query_cov"]  = len(df)
    df = df[df["span_to_scaffold_ratio"]<= args.max_span_ratio]; counts["span_ratio"] = len(df)
    df = df[df["max_mapq"]              >= args.min_mapq];        counts["mapq"]       = len(df)
    df = df[df["n_alignment_blocks"]    <= args.max_blocks];      counts["blocks"]     = len(df)

    counts["passing"] = len(df)
    return df, counts


# ── REPORT ────────────────────────────────────────────────────────────────────

def print_summary(counts, args):
    print("[filter_scaffolds] Thresholds applied:", flush=True)
    print(f"  identity_pct          >= {args.min_identity}", flush=True)
    print(f"  query_coverage_pct    >= {args.min_query_cov}", flush=True)
    print(f"  span_to_scaffold_ratio<= {args.max_span_ratio}", flush=True)
    print(f"  max_mapq              >= {args.min_mapq}", flush=True)
    print(f"  n_alignment_blocks    <= {args.max_blocks}", flush=True)
    print(f"[filter_scaffolds] Input scaffolds:      {counts['input']}", flush=True)
    print(f"[filter_scaffolds] After identity:       {counts['identity']}", flush=True)
    print(f"[filter_scaffolds] After query coverage: {counts['query_cov']}", flush=True)
    print(f"[filter_scaffolds] After span ratio:     {counts['span_ratio']}", flush=True)
    print(f"[filter_scaffolds] After MAPQ:           {counts['mapq']}", flush=True)
    print(f"[filter_scaffolds] After block count:    {counts['blocks']}", flush=True)
    print(f"[filter_scaffolds] Passing scaffolds:    {counts['passing']}", flush=True)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"[filter_scaffolds] ERROR: input not found: {args.input}")

    df = pd.read_csv(args.input, sep="\t", dtype={"best_chromosome": str})

    # dominance_ratio column contains "inf" strings — convert to float for potential future use
    # (not filtered here, but keep numeric where possible)
    df["dominance_ratio"] = pd.to_numeric(df["dominance_ratio"], errors="coerce")

    passing, counts = filter_scaffolds(df, args)
    print_summary(counts, args)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    passing.to_csv(args.output, sep="\t", index=False)
    print(f"[filter_scaffolds] Written {args.output}", flush=True)

    # Per-chromosome breakdown
    if len(passing) > 0:
        print("[filter_scaffolds] Passing scaffolds per chromosome:", flush=True)
        for chrom, n in passing.groupby("best_chromosome").size().sort_index().items():
            print(f"  {chrom}: {n}", flush=True)


if __name__ == "__main__":
    main()
