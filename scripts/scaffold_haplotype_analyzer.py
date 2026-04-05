"""
scaffold_haplotype_analyzer.py — Characterize unplaced scaffolds as alt-haplotypes or repeats.

Partitions a reference genome into placed chromosomes and unplaced scaffolds, aligns all
scaffolds to the chromosomal reference with minimap2 -x asm5, and reports per-scaffold
alignment statistics to identify candidate alternative haplotypes.

Completely independent of the remapper pipeline — takes only a reference FASTA.

Usage:
  python scripts/scaffold_haplotype_analyzer.py \\
      -r equCab3/equCab3_genome.fa \\
      -o remap_assessment/scaffold_haplotype_analysis \\
      [--scaffold-pattern Un_NW_] \\
      [--threads 8] \\
      [--keep-temp]
"""

import argparse
import math
import os
import re
import subprocess
import sys
from collections import defaultdict


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Characterize unplaced scaffolds vs. placed chromosomes via asm5 alignment."
    )
    p.add_argument("-r", "--reference",  required=True,
                   help="Reference genome FASTA (must have .fai beside it)")
    p.add_argument("-o", "--output-dir", required=True,
                   help="Output directory (created if absent)")
    p.add_argument("--scaffold-pattern", default="Un_NW_",
                   help="Substring identifying unplaced scaffold names in FAI (default: Un_NW_)")
    p.add_argument("-t", "--threads", type=int, default=4,
                   help="Threads for minimap2 (default: 4)")
    p.add_argument("--keep-temp", action="store_true",
                   help="Retain intermediate FASTA files after completion")
    return p.parse_args()


# ── FAI PARTITION ─────────────────────────────────────────────────────────────

def partition_sequences(fai_path, scaffold_pattern):
    """
    Read FAI; split sequence names into scaffold_ids and chrom_ids.
    Returns (scaffold_ids: list, chrom_ids: list).
    """
    scaffold_ids = []
    chrom_ids    = []
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            if scaffold_pattern in name:
                scaffold_ids.append(name)
            else:
                chrom_ids.append(name)
    return scaffold_ids, chrom_ids


# ── SEQUENCE EXTRACTION ───────────────────────────────────────────────────────

def extract_sequences(seq_ids, ref_fa, out_fa):
    """
    Run: samtools faidx ref_fa id1 id2 ... > out_fa
    All IDs in one call (well within Linux ARG_MAX).
    """
    cmd = ["samtools", "faidx", ref_fa] + seq_ids
    print(f"[scaffold_analysis] Extracting {len(seq_ids)} sequences to {out_fa} ...", flush=True)
    with open(out_fa, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if result.returncode != 0:
        sys.exit(
            f"[scaffold_analysis] ERROR: samtools faidx failed:\n{result.stderr.decode()}"
        )


# ── TOOL CHECKS ───────────────────────────────────────────────────────────────

def _check_tools():
    """Exit with a clear message if samtools or minimap2 are not on PATH."""
    for tool in ("samtools", "minimap2"):
        result = subprocess.run(["which", tool], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            sys.exit(f"[scaffold_analysis] ERROR: '{tool}' not found on PATH. "
                     f"Activate the 'remap' conda environment.")


# ── ALIGNMENT ─────────────────────────────────────────────────────────────────

def run_minimap2_asm5(ref_fa, query_fa, paf_path, threads):
    """
    Align query_fa against ref_fa using minimap2 asm5 preset.
    --cs  : output cs tag (for indel counting)
    -c    : output CIGAR in PAF
    """
    cmd = [
        "minimap2", "-x", "asm5",
        "--cs", "-c",
        "-t", str(threads),
        ref_fa, query_fa,
    ]
    print(f"[scaffold_analysis] Running minimap2 asm5 ({threads} threads) ...", flush=True)
    with open(paf_path, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if result.returncode != 0:
        sys.exit(
            f"[scaffold_analysis] ERROR: minimap2 failed:\n{result.stderr.decode()}"
        )
    print(f"[scaffold_analysis] PAF written to {paf_path}", flush=True)


# ── PAF HELPERS ───────────────────────────────────────────────────────────────

def _count_cs_indels(fields):
    """
    Scan optional PAF fields for the cs:Z: tag.
    Count insertion events (+[acgt]+) and deletion events (-[acgt]+).
    minimap2 asm5 cs tag uses lowercase ACGT.
    Returns (n_insertions, n_deletions).
    """
    cs_str = ""
    for f in fields:
        if f.startswith("cs:Z:"):
            cs_str = f[5:]
            break
    n_ins = len(re.findall(r"\+[acgt]+", cs_str))
    n_del = len(re.findall(r"-[acgt]+", cs_str))
    return n_ins, n_del


def _union_length(intervals):
    """Merge overlapping (start, end) pairs; return total covered bases."""
    if not intervals:
        return 0
    merged = sorted(intervals, key=lambda x: x[0])
    total = 0
    cur_start, cur_end = merged[0]
    for s, e in merged[1:]:
        if s <= cur_end:
            cur_end = max(cur_end, e)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = s, e
    total += cur_end - cur_start
    return total


# ── PAF PARSING ───────────────────────────────────────────────────────────────

def parse_paf(paf_path):
    """
    Parse PAF file. Returns dict[scaffold_id -> list[record_dict]].
    No tname filtering needed — the target is the chromosomal-only reference,
    so every tname is guaranteed to be a placed chromosome.

    PAF columns (0-based):
      0  qname   1  qlen   2  qstart  3  qend   4  strand
      5  tname   6  tlen   7  tstart  8  tend    9  nmatch
      10 alen    11 mapq   12+ optional fields
    """
    hits = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue

            qname  = parts[0]
            qstart = int(parts[2])
            qend   = int(parts[3])
            strand = parts[4]
            tname  = parts[5]
            tstart = int(parts[7])
            tend   = int(parts[8])
            nmatch = int(parts[9])
            alen   = int(parts[10])
            mapq   = int(parts[11])

            n_ins, n_del = _count_cs_indels(parts[12:])

            hits[qname].append({
                "qstart": qstart,
                "qend":   qend,
                "strand": strand,
                "tname":  tname,
                "tstart": tstart,
                "tend":   tend,
                "nmatch": nmatch,
                "alen":   alen,
                "mapq":   mapq,
                "n_ins":  n_ins,
                "n_del":  n_del,
            })
    return dict(hits)


# ── METRIC AGGREGATION ────────────────────────────────────────────────────────

def compute_scaffold_metrics(scaffold_id, scaffold_len, paf_records):
    """
    Aggregate PAF records for one scaffold into 13 output metrics.

    Best chromosome = chromosome with highest total nmatch across all its PAF records.

    Returns a dict with keys matching scaffold_summary.tsv columns, or None if
    paf_records is empty (caller should skip empty scaffolds from output).
    """
    if not paf_records:
        return None

    # ── Group records by target chromosome ───────────────────────────────────
    by_chrom = defaultdict(list)
    for rec in paf_records:
        by_chrom[rec["tname"]].append(rec)

    # ── Per-chromosome aggregates ─────────────────────────────────────────────
    chrom_stats = {}
    for chrom, recs in by_chrom.items():
        total_nm  = sum(r["nmatch"] for r in recs)
        total_al  = sum(r["alen"]   for r in recs)
        covered   = _union_length([(r["qstart"], r["qend"]) for r in recs])
        n_indels  = sum(r["n_ins"] + r["n_del"] for r in recs)
        tstart_min = min(r["tstart"] for r in recs)
        tend_max   = max(r["tend"]   for r in recs)
        max_mapq   = max(r["mapq"]   for r in recs)
        n_blocks   = len(recs)
        chrom_stats[chrom] = {
            "nmatch":      total_nm,
            "alen":        total_al,
            "covered":     covered,
            "n_indels":    n_indels,
            "tstart_min":  tstart_min,
            "tend_max":    tend_max,
            "max_mapq":    max_mapq,
            "n_blocks":    n_blocks,
        }

    # ── Rank chromosomes by total nmatch ─────────────────────────────────────
    ranked = sorted(chrom_stats.keys(), key=lambda c: chrom_stats[c]["nmatch"], reverse=True)
    best   = ranked[0]
    bs     = chrom_stats[best]

    # ── Derived metrics for best chromosome ──────────────────────────────────
    identity_pct       = bs["nmatch"] / bs["alen"] * 100 if bs["alen"] > 0 else 0.0
    query_coverage_pct = bs["covered"] / scaffold_len * 100 if scaffold_len > 0 else 0.0
    target_span_bp     = bs["tend_max"] - bs["tstart_min"]
    span_to_scaffold   = target_span_bp / scaffold_len if scaffold_len > 0 else 0.0

    if len(ranked) >= 2:
        second_nm    = chrom_stats[ranked[1]]["nmatch"]
        dominance    = bs["nmatch"] / second_nm if second_nm > 0 else math.inf
    else:
        dominance    = math.inf

    # ── all_chromosome_hits string ────────────────────────────────────────────
    parts = []
    for c in ranked:
        cov_pct = chrom_stats[c]["covered"] / scaffold_len * 100 if scaffold_len > 0 else 0.0
        parts.append(f"{c}:nmatch={chrom_stats[c]['nmatch']},cov={cov_pct:.1f}%")
    all_hits = ";".join(parts)

    return {
        "scaffold_id":          scaffold_id,
        "scaffold_len":         scaffold_len,
        "best_chromosome":      best,
        "identity_pct":         round(identity_pct, 3),
        "query_coverage_pct":   round(query_coverage_pct, 3),
        "n_indels":             bs["n_indels"],
        "dominance_ratio":      round(dominance, 3) if not math.isinf(dominance) else "inf",
        "n_chr_hits":           len(ranked),
        "n_alignment_blocks":   bs["n_blocks"],
        "target_span_bp":       target_span_bp,
        "span_to_scaffold_ratio": round(span_to_scaffold, 3),
        "max_mapq":             bs["max_mapq"],
        "all_chromosome_hits":  all_hits,
    }

# ── OUTPUT WRITERS ────────────────────────────────────────────────────────────

SUMMARY_COLS = [
    "scaffold_id", "scaffold_len", "best_chromosome",
    "identity_pct", "query_coverage_pct", "n_indels",
    "dominance_ratio", "n_chr_hits", "n_alignment_blocks",
    "target_span_bp", "span_to_scaffold_ratio", "max_mapq",
    "all_chromosome_hits",
]


def write_scaffold_summary(results, output_dir):
    """Write scaffold_summary.tsv — one row per scaffold with ≥1 chromosome hit."""
    path = os.path.join(output_dir, "scaffold_summary.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(SUMMARY_COLS) + "\n")
        for r in sorted(results, key=lambda x: x["scaffold_id"]):
            fh.write("\t".join(str(r[c]) for c in SUMMARY_COLS) + "\n")
    print(f"[scaffold_analysis] Written {path} ({len(results)} rows)", flush=True)


def write_report(n_total, results, output_dir):
    """
    Write condensed report.txt.
    n_total: total scaffolds tested (from FAI).
    results: list of metric dicts (scaffolds with >=1 hit only).
    """
    n_hits   = len(results)
    n_no_hit = n_total - n_hits

    # Identity buckets
    id_buckets = {">=99%": 0, "95-99%": 0, "90-95%": 0, "<90%": 0}
    for r in results:
        v = float(r["identity_pct"])
        if   v >= 99: id_buckets[">=99%"]   += 1
        elif v >= 95: id_buckets["95-99%"]  += 1
        elif v >= 90: id_buckets["90-95%"]  += 1
        else:         id_buckets["<90%"]    += 1

    # span_to_scaffold_ratio buckets
    span_buckets = {"<=1.5": 0, "1.5-5": 0, ">5": 0}
    for r in results:
        v = float(r["span_to_scaffold_ratio"])
        if   v <= 1.5: span_buckets["<=1.5"] += 1
        elif v <= 5.0: span_buckets["1.5-5"] += 1
        else:          span_buckets[">5"]    += 1

    # max_mapq buckets
    mapq_buckets = {"60": 0, "1-59": 0, "0": 0}
    for r in results:
        v = int(r["max_mapq"])
        if   v == 60: mapq_buckets["60"]   += 1
        elif v >= 1:  mapq_buckets["1-59"] += 1
        else:         mapq_buckets["0"]    += 1

    # Per-chromosome counts
    chrom_counts = defaultdict(int)
    for r in results:
        chrom_counts[r["best_chromosome"]] += 1

    lines = [
        "=" * 70,
        "Scaffold Haplotype Analyzer — Report",
        "=" * 70,
        "",
        f"Total scaffolds tested:          {n_total}",
        f"Scaffolds with >=1 chr hit:      {n_hits}",
        f"Scaffolds with zero chr hits:    {n_no_hit}",
        "",
        "Identity distribution (scaffolds with >=1 hit):",
    ]
    for bucket, count in id_buckets.items():
        lines.append(f"  {bucket}: {count}")
    lines += [
        "",
        "span_to_scaffold_ratio distribution:",
    ]
    for bucket, count in span_buckets.items():
        lines.append(f"  {bucket}: {count}")
    lines += [
        "",
        "max_mapq distribution:",
    ]
    for bucket, count in mapq_buckets.items():
        lines.append(f"  MAPQ={bucket}: {count}")
    lines += [
        "",
        "Best-chromosome assignment counts:",
    ]
    for chrom in sorted(chrom_counts.keys()):
        lines.append(f"  {chrom}: {chrom_counts[chrom]}")

    path = os.path.join(output_dir, "report.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[scaffold_analysis] Written {path}", flush=True)


# ── CLEANUP ───────────────────────────────────────────────────────────────────

def cleanup_temp(paths):
    """Remove temp files."""
    for p in paths:
        if os.path.exists(p):
            os.remove(p)
            print(f"[scaffold_analysis] Removed temp file: {p}", flush=True)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── Validate inputs ───────────────────────────────────────────────────────
    fai_path = args.reference + ".fai"
    if not os.path.exists(args.reference):
        sys.exit(f"[scaffold_analysis] ERROR: reference not found: {args.reference}")
    if not os.path.exists(fai_path):
        sys.exit(f"[scaffold_analysis] ERROR: FAI not found: {fai_path}  "
                 f"Run: samtools faidx {args.reference}")
    _check_tools()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Partition FAI ─────────────────────────────────────────────────────────
    print(f"[scaffold_analysis] Partitioning FAI with pattern '{args.scaffold_pattern}' ...",
          flush=True)
    scaffold_ids, chrom_ids = partition_sequences(fai_path, args.scaffold_pattern)

    if not scaffold_ids:
        sys.exit(f"[scaffold_analysis] ERROR: no sequences match "
                 f"--scaffold-pattern '{args.scaffold_pattern}'. "
                 f"Check that the pattern is correct for this reference.")
    if not chrom_ids:
        sys.exit(f"[scaffold_analysis] ERROR: no chromosome sequences remain after "
                 f"partitioning. The reference appears to contain only scaffold sequences.")

    print(f"[scaffold_analysis] {len(scaffold_ids)} scaffolds, "
          f"{len(chrom_ids)} chromosomes", flush=True)

    # ── Extract FASTAs ────────────────────────────────────────────────────────
    scaffold_fa = os.path.join(args.output_dir, "unplaced_scaffolds.fa")
    chrom_fa    = os.path.join(args.output_dir, "chromosomal_ref.fa")
    extract_sequences(scaffold_ids, args.reference, scaffold_fa)
    extract_sequences(chrom_ids,    args.reference, chrom_fa)

    # ── Align ─────────────────────────────────────────────────────────────────
    paf_path = os.path.join(args.output_dir, "scaffold_vs_chr.paf")
    run_minimap2_asm5(chrom_fa, scaffold_fa, paf_path, args.threads)

    # ── Parse PAF ─────────────────────────────────────────────────────────────
    print("[scaffold_analysis] Parsing PAF ...", flush=True)
    paf_hits = parse_paf(paf_path)
    n_with_hits = len(paf_hits)
    print(f"[scaffold_analysis] {n_with_hits}/{len(scaffold_ids)} scaffolds have chromosome hits",
          flush=True)

    # ── Compute metrics ───────────────────────────────────────────────────────
    print("[scaffold_analysis] Computing per-scaffold metrics ...", flush=True)

    # Build scaffold_len lookup from FAI
    scaffold_len = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 2 and args.scaffold_pattern in parts[0]:
                scaffold_len[parts[0]] = int(parts[1])

    results = []
    for sid in scaffold_ids:
        recs    = paf_hits.get(sid, [])
        metrics = compute_scaffold_metrics(sid, scaffold_len.get(sid, 0), recs)
        if metrics is not None:
            results.append(metrics)

    # ── Write outputs ─────────────────────────────────────────────────────────
    print("[scaffold_analysis] Writing outputs ...", flush=True)
    write_scaffold_summary(results, args.output_dir)
    write_report(len(scaffold_ids), results, args.output_dir)

    # ── Cleanup ───────────────────────────────────────────────────────────────
    if not args.keep_temp:
        cleanup_temp([scaffold_fa, chrom_fa])

    print(f"[scaffold_analysis] Done. Outputs in: {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
