"""
analyze_unplaced_scaffolds.py — Characterize unplaced scaffold markers in remapped manifests.

For markers landing on Un_NW_* scaffolds, aligns each scaffold to the full reference with
minimap2 -x asm5 to determine whether the scaffold is an alternative haplotype of a placed
chromosome or an independent duplication/repeat.

Usage:
  python scripts/analyze_unplaced_scaffolds.py \\
      -i results/Equine80select_24_20067593_B1_remapped_equCab3.csv \\
      -r equCab3/equCab3_genome.fa \\
      -o remap_assessment/unplaced_scaffold_analysis \\
      --assembly equCab3 \\
      --threads 8
"""

import argparse
import math
import os
import re
import subprocess
import sys
from collections import defaultdict

import pandas as pd

try:
    from scipy import stats as scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Characterize unplaced scaffold markers vs. placed chromosomes."
    )
    p.add_argument("-i", "--input",      required=True, help="Remapped manifest CSV")
    p.add_argument("-r", "--reference",  required=True, help="Reference genome FASTA")
    p.add_argument("-o", "--output-dir", required=True, help="Output directory")
    p.add_argument(
        "-a", "--assembly",
        default="equCab3",
        help="Assembly name used in CSV column names (default: equCab3)",
    )
    p.add_argument("-t", "--threads", type=int, default=4, help="Threads for minimap2 (default: 4)")
    return p.parse_args()


# ── STEP 1: PARSE MARKERS ────────────────────────────────────────────────────

def parse_markers_from_csv(csv_path, assembly):
    """
    Load the remapped manifest and return:
      - scaffold_markers: dict[scaffold_id → list[marker_dicts]]
      - all scaffolds that contain 'Un_NW_'
    """
    col_chr   = f"Chr_{assembly}"
    col_pos   = f"MapInfo_{assembly}"
    col_mapq_ts  = "MAPQ_TopGenomicSeq"
    col_mapq_p   = "MAPQ_Probe"

    df = pd.read_csv(csv_path, dtype={col_chr: str}, low_memory=False)

    mask = df[col_chr].str.contains("Un_NW_", na=False)
    unplaced = df[mask].copy()

    scaffold_markers = defaultdict(list)
    for _, row in unplaced.iterrows():
        scaffold_id = row[col_chr]
        marker = {
            "marker_id":        row.get("Name", row.get("IlmnID", "unknown")),
            "scaffold_id":      scaffold_id,
            "scaffold_pos":     int(row[col_pos]),
            "mapq_topgenomic":  int(row[col_mapq_ts]),
            "mapq_probe":       int(row[col_mapq_p]),
        }
        scaffold_markers[scaffold_id].append(marker)

    return dict(scaffold_markers)


# ── STEP 2: READ FAI ─────────────────────────────────────────────────────────

def read_fai_sizes(fai_path):
    """Parse FAI file; return dict[seq_name → length]."""
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


# ── STEP 3: EXTRACT SCAFFOLD SEQUENCES ───────────────────────────────────────

def extract_scaffold_sequences(scaffold_ids, ref_fa, out_fa):
    """
    Run: samtools faidx ref.fa scaf1 scaf2 ... > out_fa
    All scaffold IDs in one call (well within Linux ARG_MAX).
    """
    cmd = ["samtools", "faidx", ref_fa] + sorted(scaffold_ids)
    print(f"[scaffold_analysis] Extracting {len(scaffold_ids)} scaffold sequences ...", flush=True)
    with open(out_fa, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if result.returncode != 0:
        sys.exit(f"[scaffold_analysis] ERROR: samtools faidx failed:\n{result.stderr.decode()}")
    print(f"[scaffold_analysis] Scaffold FASTA written to {out_fa}", flush=True)


# ── STEP 4: RUN MINIMAP2 ASM5 ────────────────────────────────────────────────

def run_minimap2_asm5(ref_fa, query_fa, paf_path, threads):
    """
    Align scaffold sequences to the full reference using asm5 preset.
    --cs tag for cs string (indel counting).
    -c for CIGAR in PAF output.
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
        sys.exit(f"[scaffold_analysis] ERROR: minimap2 failed:\n{result.stderr.decode()}")
    print(f"[scaffold_analysis] PAF written to {paf_path}", flush=True)


# ── STEP 5: PARSE PAF ────────────────────────────────────────────────────────

def _count_cs_indels(fields):
    """
    Scan optional PAF fields for the cs:Z: tag.
    Count insertion events (+[acgt]+) and deletion events (-[acgt]+).
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


def parse_paf(paf_path):
    """
    Parse PAF file. Skip records where tname contains 'Un_NW_' (scaffold-vs-scaffold).
    Returns dict[query_name → list[PafRecord]] where PafRecord is a dict.
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
            tname = parts[5]
            if "Un_NW_" in tname:
                continue

            qname    = parts[0]
            qlen     = int(parts[1])
            qstart   = int(parts[2])
            qend     = int(parts[3])
            strand   = parts[4]
            tlen     = int(parts[6])
            tstart   = int(parts[7])
            tend     = int(parts[8])
            nmatch   = int(parts[9])
            alen     = int(parts[10])
            mapq     = int(parts[11])

            n_ins, n_del = _count_cs_indels(parts[12:])

            hits[qname].append({
                "qname":   qname,
                "qlen":    qlen,
                "qstart":  qstart,
                "qend":    qend,
                "strand":  strand,
                "tname":   tname,
                "tlen":    tlen,
                "tstart":  tstart,
                "tend":    tend,
                "nmatch":  nmatch,
                "alen":    alen,
                "mapq":    mapq,
                "n_ins":   n_ins,
                "n_del":   n_del,
            })
    return dict(hits)


# ── STEP 6: COMPUTE METRICS ──────────────────────────────────────────────────

def _union_length(intervals):
    """Merge overlapping (start, end) tuples; return total covered length."""
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


def project_marker_to_chromosome(marker_pos, paf_records_for_best_chrom):
    """
    Find the PAF block on the best chromosome that covers (marker_pos - 1) in query coordinates.
    Returns (chrom_pos, block) or (None, None).

    marker_pos is 1-based; PAF qstart/qend are 0-based half-open.
    """
    marker_0based = marker_pos - 1
    for rec in paf_records_for_best_chrom:
        if rec["qstart"] <= marker_0based < rec["qend"]:
            if rec["strand"] == "+":
                chrom_pos = rec["tstart"] + (marker_0based - rec["qstart"]) + 1
            else:
                chrom_pos = rec["tend"] - (marker_0based - rec["qstart"])
            return chrom_pos, rec
    return None, None


def compute_scaffold_metrics(scaffold_id, scaffold_len, paf_records, markers):
    """
    Aggregate per-chromosome metrics for one scaffold.
    Returns a dict of metrics.
    """
    # Group PAF records by target chromosome
    by_chrom = defaultdict(list)
    for rec in paf_records:
        by_chrom[rec["tname"]].append(rec)

    if not by_chrom:
        return {
            "scaffold_id":         scaffold_id,
            "scaffold_len":        scaffold_len,
            "n_markers":           len(markers),
            "best_chromosome":     "NA",
            "identity_pct":        "NA",
            "query_coverage_pct":  "NA",
            "n_indels":            "NA",
            "dominance_ratio":     "NA",
            "n_chr_hits":          0,
            "synteny_corr":        "NA",
            "synteny_pval":        "NA",
            "all_chromosome_hits": "NA",
        }

    # Aggregate nmatch per chromosome (used for dominance and ordering)
    chrom_nmatch = {}
    chrom_alen   = {}
    chrom_qcov   = {}
    chrom_indels = {}

    for chrom, recs in by_chrom.items():
        total_nm = sum(r["nmatch"] for r in recs)
        total_al = sum(r["alen"]   for r in recs)
        intervals = [(r["qstart"], r["qend"]) for r in recs]
        covered   = _union_length(intervals)
        n_indels  = sum(r["n_ins"] + r["n_del"] for r in recs)

        chrom_nmatch[chrom] = total_nm
        chrom_alen[chrom]   = total_al
        chrom_qcov[chrom]   = covered
        chrom_indels[chrom] = n_indels

    # Sort chromosomes by total nmatch descending
    ranked = sorted(chrom_nmatch.keys(), key=lambda c: chrom_nmatch[c], reverse=True)
    best_chrom  = ranked[0]
    best_nmatch = chrom_nmatch[best_chrom]
    best_alen   = chrom_alen[best_chrom]
    best_qcov   = chrom_qcov[best_chrom]
    best_indels = chrom_indels[best_chrom]

    identity_pct       = best_nmatch / best_alen * 100  if best_alen > 0 else 0.0
    query_coverage_pct = best_qcov   / scaffold_len * 100 if scaffold_len > 0 else 0.0

    # Dominance ratio
    if len(ranked) >= 2:
        second_nmatch = chrom_nmatch[ranked[1]]
        dominance_ratio = best_nmatch / second_nmatch if second_nmatch > 0 else math.inf
    else:
        dominance_ratio = math.inf

    # Build all_chromosome_hits string (sorted by nmatch desc)
    all_hits_parts = []
    for c in ranked:
        cov_pct = chrom_qcov[c] / scaffold_len * 100 if scaffold_len > 0 else 0.0
        all_hits_parts.append(f"{c}:nmatch={chrom_nmatch[c]},cov={cov_pct:.1f}%")
    all_chromosome_hits = ";".join(all_hits_parts)

    # Synteny check (requires ≥2 markers and ≥1 PAF hit on best chrom)
    best_recs = by_chrom[best_chrom]
    synteny_corr = "NA"
    synteny_pval = "NA"

    if len(markers) >= 2 and best_recs:
        scaffold_positions = []
        chrom_positions    = []
        for m in markers:
            cp, _ = project_marker_to_chromosome(m["scaffold_pos"], best_recs)
            if cp is not None:
                scaffold_positions.append(m["scaffold_pos"])
                chrom_positions.append(cp)

        if len(scaffold_positions) >= 2:
            if HAS_SCIPY:
                rho, pval = scipy_stats.spearmanr(scaffold_positions, chrom_positions)
                synteny_corr = f"{rho:.4f}"
                synteny_pval = f"{pval:.4g}" if not math.isnan(pval) else "nan"
            else:
                # Manual rank correlation if scipy unavailable
                n = len(scaffold_positions)
                def ranks(lst):
                    sorted_lst = sorted(enumerate(lst), key=lambda x: x[1])
                    r = [0] * n
                    for rank, (idx, _) in enumerate(sorted_lst):
                        r[idx] = rank + 1
                    return r
                rs = ranks(scaffold_positions)
                rc = ranks(chrom_positions)
                d2 = sum((a - b) ** 2 for a, b in zip(rs, rc))
                rho = 1 - 6 * d2 / (n * (n * n - 1))
                synteny_corr = f"{rho:.4f}"
                synteny_pval = "NA(no_scipy)"

    return {
        "scaffold_id":         scaffold_id,
        "scaffold_len":        scaffold_len,
        "n_markers":           len(markers),
        "best_chromosome":     best_chrom,
        "identity_pct":        f"{identity_pct:.3f}",
        "query_coverage_pct":  f"{query_coverage_pct:.3f}",
        "n_indels":            best_indels,
        "dominance_ratio":     f"{dominance_ratio:.3f}" if not math.isinf(dominance_ratio) else "inf",
        "n_chr_hits":          len(ranked),
        "synteny_corr":        synteny_corr,
        "synteny_pval":        synteny_pval,
        "all_chromosome_hits": all_chromosome_hits,
    }


# ── STEP 7: WRITE OUTPUTS ────────────────────────────────────────────────────

def write_outputs(results, scaffold_markers, paf_hits, output_dir):
    """Write scaffold_summary.tsv, marker_summary.tsv, and report.txt."""

    # ── scaffold_summary.tsv ─────────────────────────────────────────────────
    scaffold_tsv = os.path.join(output_dir, "scaffold_summary.tsv")
    scaffold_cols = [
        "scaffold_id", "scaffold_len", "n_markers",
        "best_chromosome", "identity_pct", "query_coverage_pct",
        "n_indels", "dominance_ratio", "n_chr_hits",
        "synteny_corr", "synteny_pval", "all_chromosome_hits",
    ]
    with open(scaffold_tsv, "w") as fh:
        fh.write("\t".join(scaffold_cols) + "\n")
        for r in sorted(results, key=lambda x: x["scaffold_id"]):
            row = [str(r[c]) for c in scaffold_cols]
            fh.write("\t".join(row) + "\n")
    print(f"[scaffold_analysis] Written {scaffold_tsv}", flush=True)

    # ── marker_summary.tsv ───────────────────────────────────────────────────
    marker_tsv = os.path.join(output_dir, "marker_summary.tsv")
    marker_cols = [
        "marker_id", "scaffold_id", "scaffold_pos",
        "mapq_topgenomic", "mapq_probe",
        "best_chromosome", "projected_chrom_pos",
        "covering_block_qstart", "covering_block_qend",
        "covering_block_tstart", "covering_block_tend", "covering_block_strand",
    ]
    # Build best_chrom lookup from results
    best_chrom_for = {r["scaffold_id"]: r["best_chromosome"] for r in results}

    with open(marker_tsv, "w") as fh:
        fh.write("\t".join(marker_cols) + "\n")
        for scaffold_id, markers in sorted(scaffold_markers.items()):
            best_chrom = best_chrom_for.get(scaffold_id, "NA")
            recs_for_best = []
            if best_chrom != "NA":
                recs_for_best = [r for r in paf_hits.get(scaffold_id, [])
                                 if r["tname"] == best_chrom]
            for m in sorted(markers, key=lambda x: x["scaffold_pos"]):
                cp, block = project_marker_to_chromosome(m["scaffold_pos"], recs_for_best)
                row = [
                    m["marker_id"],
                    scaffold_id,
                    str(m["scaffold_pos"]),
                    str(m["mapq_topgenomic"]),
                    str(m["mapq_probe"]),
                    best_chrom,
                    str(cp) if cp is not None else "NA",
                    str(block["qstart"])  if block else "NA",
                    str(block["qend"])    if block else "NA",
                    str(block["tstart"])  if block else "NA",
                    str(block["tend"])    if block else "NA",
                    block["strand"]       if block else "NA",
                ]
                fh.write("\t".join(row) + "\n")
    print(f"[scaffold_analysis] Written {marker_tsv}", flush=True)

    # ── report.txt ───────────────────────────────────────────────────────────
    report_path = os.path.join(output_dir, "report.txt")
    _write_report(results, scaffold_markers, report_path)
    print(f"[scaffold_analysis] Written {report_path}", flush=True)


def _write_report(results, scaffold_markers, report_path):
    total_scaffolds = len(results)
    total_markers   = sum(r["n_markers"] for r in results)

    no_hit_scaffolds = [r for r in results if r["n_chr_hits"] == 0]

    # MAPQ distribution (from markers)
    mapq_dist = defaultdict(int)
    for markers in scaffold_markers.values():
        for m in markers:
            mapq_dist[m["mapq_topgenomic"]] += 1

    # scaffolds with chromosome hits
    has_hit = [r for r in results if r["n_chr_hits"] > 0]

    # Per-chromosome scaffold counts (best chromosome)
    chrom_counts = defaultdict(int)
    for r in has_hit:
        chrom_counts[r["best_chromosome"]] += 1

    # Identity distribution buckets (for scaffolds with hits)
    id_buckets = {">=99%": 0, "95-99%": 0, "90-95%": 0, "<90%": 0}
    for r in has_hit:
        try:
            pct = float(r["identity_pct"])
        except (ValueError, TypeError):
            continue
        if pct >= 99:
            id_buckets[">=99%"] += 1
        elif pct >= 95:
            id_buckets["95-99%"] += 1
        elif pct >= 90:
            id_buckets["90-95%"] += 1
        else:
            id_buckets["<90%"] += 1

    # Synteny distribution
    syn_buckets = {">=0.9": 0, "0.5-0.9": 0, "<0.5": 0, "NA": 0}
    for r in results:
        sc = r["synteny_corr"]
        if sc == "NA" or "NA" in str(sc):
            syn_buckets["NA"] += 1
        else:
            try:
                v = float(sc)
                if v >= 0.9:
                    syn_buckets[">=0.9"] += 1
                elif v >= 0.5:
                    syn_buckets["0.5-0.9"] += 1
                else:
                    syn_buckets["<0.5"] += 1
            except ValueError:
                syn_buckets["NA"] += 1

    lines = []
    lines.append("=" * 70)
    lines.append("Unplaced Scaffold Analysis — Report")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"Total unplaced scaffolds with markers: {total_scaffolds}")
    lines.append(f"Total unplaced markers:                {total_markers}")
    lines.append(f"Scaffolds with no chromosome hit:      {len(no_hit_scaffolds)}")
    lines.append(f"Scaffolds with >=1 chromosome hit:     {len(has_hit)}")
    lines.append("")

    lines.append("MAPQ_TopGenomicSeq distribution (unplaced markers):")
    for q in sorted(mapq_dist.keys()):
        lines.append(f"  MAPQ={q:3d}: {mapq_dist[q]:5d} markers")
    lines.append("")

    lines.append("Identity distribution (scaffolds with >=1 chr hit):")
    for bucket, count in id_buckets.items():
        lines.append(f"  {bucket}: {count}")
    lines.append("")

    lines.append("Synteny correlation distribution:")
    for bucket, count in syn_buckets.items():
        lines.append(f"  rho {bucket}: {count} scaffolds")
    lines.append("")

    lines.append("Best-chromosome hit counts (scaffolds assigned to each chromosome):")
    for chrom in sorted(chrom_counts.keys()):
        lines.append(f"  {chrom}: {chrom_counts[chrom]} scaffold(s)")
    lines.append("")

    if len(no_hit_scaffolds) > 0:
        lines.append(f"No-hit scaffolds ({len(no_hit_scaffolds)}):")
        for r in sorted(no_hit_scaffolds, key=lambda x: x["scaffold_id"]):
            lines.append(f"  {r['scaffold_id']} (len={r['scaffold_len']}, markers={r['n_markers']})")
        lines.append("")

    if len(no_hit_scaffolds) / total_scaffolds > 0.9 if total_scaffolds > 0 else False:
        lines.append("WARNING: >90% of scaffolds have no chromosome hit.")
        lines.append("         This may indicate a reference mismatch or alignment issue.")
        lines.append("")

    with open(report_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Parse markers
    print("[scaffold_analysis] Parsing unplaced markers from CSV ...", flush=True)
    scaffold_markers = parse_markers_from_csv(args.input, args.assembly)
    n_scaffolds = len(scaffold_markers)
    n_markers   = sum(len(v) for v in scaffold_markers.values())
    print(f"[scaffold_analysis] Found {n_markers} markers on {n_scaffolds} scaffolds", flush=True)

    # Step 2: Read FAI sizes
    fai_path = args.reference + ".fai"
    if not os.path.exists(fai_path):
        sys.exit(f"[scaffold_analysis] ERROR: FAI not found: {fai_path}")
    fai_sizes = read_fai_sizes(fai_path)

    missing = [s for s in scaffold_markers if s not in fai_sizes]
    if missing:
        print(f"[scaffold_analysis] WARNING: {len(missing)} scaffold(s) not found in FAI", flush=True)

    # Step 3: Extract scaffold sequences
    scaffold_fa  = os.path.join(args.output_dir, "unplaced_scaffolds.fa")
    extract_scaffold_sequences(list(scaffold_markers.keys()), args.reference, scaffold_fa)

    # Step 4: Run minimap2 asm5
    paf_path = os.path.join(args.output_dir, "scaffold_vs_chr.paf")
    run_minimap2_asm5(args.reference, scaffold_fa, paf_path, args.threads)

    # Step 5: Parse PAF
    print("[scaffold_analysis] Parsing PAF ...", flush=True)
    paf_hits = parse_paf(paf_path)
    n_with_hits = sum(1 for s in scaffold_markers if s in paf_hits)
    print(f"[scaffold_analysis] {n_with_hits}/{n_scaffolds} scaffolds have chromosome hits", flush=True)

    # Step 6: Compute metrics for each scaffold
    print("[scaffold_analysis] Computing per-scaffold metrics ...", flush=True)
    results = []
    for scaffold_id, markers in scaffold_markers.items():
        scaffold_len = fai_sizes.get(scaffold_id, 0)
        recs = paf_hits.get(scaffold_id, [])
        metrics = compute_scaffold_metrics(scaffold_id, scaffold_len, recs, markers)
        results.append(metrics)

    # Step 7: Write outputs
    print("[scaffold_analysis] Writing outputs ...", flush=True)
    write_outputs(results, scaffold_markers, paf_hits, args.output_dir)

    # Verification summary
    print(f"[scaffold_analysis] scaffold_summary.tsv rows: {len(results)} (expected: {n_scaffolds})", flush=True)
    print(f"[scaffold_analysis] marker_summary.tsv rows:   {n_markers} (expected: {n_markers})", flush=True)
    print(f"[scaffold_analysis] Done. Outputs in: {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
