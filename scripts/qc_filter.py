"""
qc_filter.py — QC filtering, allele decision, VCF/BIM/map generation.

Takes the remapped manifest (output of remap_manifest.py) and applies a cascade of
quality filters, then generates all downstream output files:

  Filter cascade (with marker counts at each stage written to QC_Report.txt):
    1. Strand=N/A filter   — remove markers where Strand_{assembly} == 'N/A' (unmapped + ambiguous)
    2. MAPQ filter         — remove markers below MAPQ thresholds
    3. Design conflict     — keep only markers where remapped Ref matches genome Ref
    4. Polymorphic sites   — remove positions with conflicting Ref/Alt assignments
    5. Consistency filter  — remove markers with inconsistent probe/topseq alignments

  Outputs:
    _matchingSNPs.vcf
    _matchingSNPs_binary.vcf
    _matchingSNPs_binary_consistantMapping.vcf
    {prefix}_remapped_{assembly}.bim
    matchingSNPs_binary_consistantMapping.{assembly}_map
    QC_Report.txt
    remap_assessment/   (MAPQ histograms, benchmark against known-assembly markers)

Usage:
  python scripts/qc_filter.py \\
      -i remapped.csv \\
      -r reference.fa \\
      -v vcf_contigs.txt \\
      -a equCab3 \\
      -o output_dir/ \\
      [--mapq-topseq 30] \\
      [--mapq-probe 0] \\
      [--temp-dir /tmp/remap] \\
      [--prefix Equine80select]
"""

import argparse
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd
import pysam

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(seq):
    return seq.translate(COMPLEMENT)


# ── CLI ──────────────────────────────────────────────────────────────────────

def _mapq_int(value):
    """argparse type that validates MAPQ is an integer in [0, 60]."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r} is not an integer.")
    if not (0 <= ivalue <= 60):
        raise argparse.ArgumentTypeError(
            f"MAPQ value {ivalue} is out of range [0, 60]."
        )
    return ivalue


def parse_args():
    p = argparse.ArgumentParser(description="QC filter and output generation for remapped manifests.")
    p.add_argument("-i", "--input",    required=True, help="Remapped manifest CSV (from remap_manifest.py)")
    p.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p.add_argument("-v", "--vcf-contigs", required=True, help="VCF contig header file")
    p.add_argument("-a", "--assembly", default="new_assembly", help="Assembly name (must match remap_manifest.py -a)")
    p.add_argument("-o", "--output-dir", default=".", help="Output directory (default: current directory)")
    p.add_argument("--mapq-topseq", type=_mapq_int, default=30,
                   help="Minimum MAPQ for TopGenomicSeq alignments (0–60, default: 30); probe_only exempt")
    p.add_argument("--mapq-probe", type=_mapq_int, default=0,
                   help="Minimum MAPQ for probe alignments (0–60, default: 0 = disabled); topseq_only exempt")
    p.add_argument("--temp-dir", default=None,
                   help="Directory for intermediate files (default: output-dir)")
    p.add_argument("--prefix", default=None,
                   help="Output file prefix (default: derived from input filename)")
    p.add_argument("--coord-delta", type=int, default=-1,
                   help="Maximum allowed CoordDelta (|probe_coord - CIGAR_coord|). "
                        "-1 = disabled (default). Any value >= 0 removes markers with "
                        "CoordDelta > threshold. topseq_only and probe_only markers "
                        "(CoordDelta=-1) pass through.")
    p.add_argument(
        "--coordinate-role", choices=["High", "Moderate", "Low"], default="Moderate",
        help="Minimum coordinate role: High=topseq_n_probe only; "
             "Moderate=also allows topseq_only (default); Low=also allows probe_only."
    )
    p.add_argument(
        "--tie-label", choices=["unique", "resolved", "avoid_scaffolds"], default="resolved",
        help="Minimum tie resolution: unique; resolved=unique+*_resolved (default); "
             "avoid_scaffolds=resolved+scaffold_resolved."
    )
    p.add_argument(
        "--refalt-conf", choices=["High", "Moderate", "Low"], default="Moderate",
        help="Minimum RefAlt confidence: High=NM_match+NM_validated; "
             "Moderate=+NM_N/A+NM_tied (default); Low=+NM_only+NM_unmatch+NM_corrected."
    )
    p.add_argument(
        "--keep-indels", action="store_true", default=False,
        help="Keep indel markers in outputs (default: False = indels excluded)."
    )
    p.add_argument(
        "--keep-polymorphic", action="store_true", default=False,
        help="Keep markers at polymorphic positions (default: False = removed)."
    )
    return p.parse_args()



# ── VCF GENERATION ───────────────────────────────────────────────────────────

def build_pos_vcf(df, vcf_contigs_path, col_chr, col_pos, pos_vcf_path):
    """Writes a VCF with N as REF/ALT for bcftools to fill in real reference alleles."""
    with open(pos_vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.3\n")
        with open(vcf_contigs_path) as vc:
            f.write(vc.read())
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in df.iterrows():
            f.write(f"{row[col_chr]}\t{row[col_pos]}\t{row['Name']}\tN\t.\t.\t.\t.\n")


def extract_ref_alleles(pos_vcf, ref_fasta, ref_vcf):
    """
    Runs bcftools norm to pull real reference alleles from the genome FASTA.
    Returns a dict {snp_name: ref_allele}.
    """
    subprocess.check_call(
        f"bcftools norm -c ws -f {ref_fasta} {pos_vcf} > {ref_vcf} 2>/dev/null",
        shell=True,
    )
    # bcftools reorders records; rebuild in positional order from pos_vcf SNP order,
    # then align by ID.
    ref_map = {}
    with open(ref_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            ref_map[cols[2]] = cols[3].upper()
    return ref_map


def strand_normalize(allele, strand):
    """Returns the allele on the + (forward) strand."""
    if strand == "-":
        return allele.translate(COMPLEMENT)[::-1]
    return allele


# ── FINAL MAP FILE ───────────────────────────────────────────────────────────

def build_final_map(df_final, assembly, map_path):
    """
    Writes matchingSNPs_binary_consistantMapping.{assembly}_map with columns:
      chr  pos  snpID  SNP_alleles  genomic_alleles  SNP_ref_allele  genomic_ref_allele  decision

    SNP_alleles are from the manifest SNP column (e.g. A,G).
    genomic_alleles are the + strand remapped alleles.
    decision ('as_is' or 'complement') is inferred by matching SNP alleles to genomic alleles —
    direct match first, then complement match. This replaces the old XOR-based approach.
    """
    col_chr    = f"Chr_{assembly}"
    col_pos    = f"MapInfo_{assembly}"
    col_strand = f"Strand_{assembly}"
    col_ref    = f"Ref_{assembly}"
    col_alt    = f"Alt_{assembly}"

    errors = 0
    lines = []

    for _, row in df_final.iterrows():
        name = row["Name"]

        # Parse SNP alleles from manifest [A/G] format
        m = re.search(r"\[(.+?)/(.+?)\]", row.get("SNP", ""))
        if not m:
            continue
        snp_a, snp_b = m.group(1), m.group(2)

        # Genomic alleles on + strand
        strand = row[col_strand]
        raw_ref = row[col_ref] if pd.notna(row[col_ref]) else ""
        raw_alt = row[col_alt] if pd.notna(row[col_alt]) else ""
        gref = strand_normalize(raw_ref, strand)
        galt = strand_normalize(raw_alt, strand)

        # Complement of SNP alleles
        snp_a_comp = complement(snp_a)
        snp_b_comp = complement(snp_b)

        chr_  = row[col_chr]
        pos   = row[col_pos]
        is_indel = bool(re.search(r"\[D/I\]|\[I/D\]", row.get("SNP", "") or ""))

        # Infer decision by matching: try as_is first, then complement
        snp_ref = None
        decision = None
        snp_alleles = None
        geno_alleles = None

        if is_indel and {snp_a, snp_b} == {"D", "I"}:
            # D/I indel markers: D = absent/shorter allele, I = present/longer allele.
            # For deletions (gref=long, galt=""): I→gref, D→galt; snp_ref=I (ref has sequence).
            # For insertions (gref="", galt=long): D→gref, I→galt; snp_ref=D (ref lacks insertion).
            if len(gref) >= len(galt):
                d_geno, i_geno = galt, gref
                snp_ref = "I"
            else:
                d_geno, i_geno = gref, galt
                snp_ref = "D"
            snp_alleles = f"{snp_a},{snp_b}"
            geno_alleles = f"{d_geno},{i_geno}" if snp_a == "D" else f"{i_geno},{d_geno}"
            decision = "indel_as_is"
        elif snp_a == gref and snp_b == galt:
            decision, snp_ref = "as_is", snp_a
            snp_alleles = f"{snp_a},{snp_b}"
            geno_alleles = f"{gref},{galt}"
        elif snp_b == gref and snp_a == galt:
            decision, snp_ref = "as_is", snp_b
            snp_alleles = f"{snp_a},{snp_b}"
            geno_alleles = f"{galt},{gref}"
        elif snp_a_comp == gref and snp_b_comp == galt:
            decision, snp_ref = "complement", snp_a
            snp_alleles = f"{snp_a},{snp_b}"
            geno_alleles = f"{gref},{galt}"
        elif snp_b_comp == gref and snp_a_comp == galt:
            decision, snp_ref = "complement", snp_b
            snp_alleles = f"{snp_a},{snp_b}"
            geno_alleles = f"{galt},{gref}"
        else:
            errors += 1
            lines.append(f"Error\t{name}\tno_match\t{snp_a},{snp_b}\t{gref},{galt}")
            continue

        if is_indel and not decision.startswith("indel_"):
            decision = "indel_" + decision

        lines.append(
            f"{chr_}\t{pos}\t{name}\t{snp_alleles}\t{geno_alleles}\t{snp_ref}\t{gref}\t{decision}"
        )

    with open(map_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return errors


# ── MAPQ HISTOGRAMS ──────────────────────────────────────────────────────────

def write_mapq_histo(values, bin_size, path):
    """Writes a histogram of MAPQ scores in [low, high, count] format."""
    if values.empty:
        return
    bins = defaultdict(int)
    for v in values:
        b = int(v // bin_size)
        bins[b] += 1
    bmin, bmax = min(bins), max(bins)
    with open(path, "w") as f:
        for i in range(bmin, bmax + 1):
            f.write(f"{i * bin_size}\t{(i + 1) * bin_size}\t{bins.get(i, 0)}\n")


# ── BENCHMARK ────────────────────────────────────────────────────────────────

def benchmark_known_assembly(df, assembly, assessment_dir):
    """
    For markers where the original assembly matches the target assembly (GenomeBuild == assembly
    major version number, e.g. '3' for EquCab3), compare original Chr:Pos to remapped Chr:Pos.
    Writes a mismatch report to remap_assessment/equCab3.mismatches.
    """
    col_chr = f"Chr_{assembly}"
    col_pos = f"MapInfo_{assembly}"

    # Extract version number from assembly name (e.g. 'equCab3' → '3')
    ver_match = re.search(r"(\d+)$", assembly)
    if not ver_match:
        return
    ver = ver_match.group(1) + ".0"

    known = df[df["GenomeBuild"].astype(str) == ver].copy()
    if known.empty:
        return

    mismatches = known[
        known["Chr"].astype(str) + ":" + known["MapInfo"].astype(str) !=
        known[col_chr].astype(str) + ":" + known[col_pos].astype(str)
    ]
    mismatch_path = os.path.join(assessment_dir, f"{assembly}.mismatches")
    mismatches[["Name", "SNP", "Chr", "MapInfo", col_chr, col_pos, f"Strand_{assembly}"]].to_csv(
        mismatch_path, sep=",", index=False, header=False
    )
    print(f"[qc] Known-assembly benchmark: {len(known):,} markers, {len(mismatches):,} mismatches → {mismatch_path}")


# ── PROBE MAPQ FILTER ────────────────────────────────────────────────────────

def apply_probe_mapq_filter(df, threshold):
    """Return df with rows removed where MAPQ_Probe is defined and below threshold.

    NaN MAPQ_Probe means no probe alignment was used (topseq_only markers) —
    these are exempt from the filter regardless of threshold.
    threshold=0 disables the filter entirely and returns df unchanged.
    """
    if threshold <= 0:
        return df
    probe_mapq = df["MAPQ_Probe"]
    probe_fail = probe_mapq.notna() & (probe_mapq < threshold)
    return df[~probe_fail]



# ── INDEL DESIGN CONFLICT CHECK ──────────────────────────────────────────────

def check_deletion_ref_match(fasta, chrom, mapinfo, gref):
    """Check whether the reference genome sequence at mapinfo matches gref.

    Used to validate deletion-ref alleles: fetches len(gref) bases from the
    reference at mapinfo (1-based) and compares to gref (case-insensitive).

    For insertion-ref alleles (gref == ''), there is nothing to verify —
    returns True (no conflict detected).

    fasta   : open pysam.FastaFile
    chrom   : chromosome name
    mapinfo : 1-based start position of the deletion sequence
    gref    : strand-normalised ref allele (empty string for insertions)

    Returns True (no conflict) or False (mismatch or fetch error).
    """
    if gref == "":
        return True  # insertion: ref is empty, nothing to verify
    try:
        ref_seq = fasta.fetch(chrom, mapinfo - 1, mapinfo - 1 + len(gref)).upper()
    except (ValueError, KeyError):
        return False
    return ref_seq == gref.upper()


# ── ANCHOR-BASE VCF ENCODING FOR INDELS ──────────────────────────────────────

def make_anchor_alleles(fasta, chrom, mapinfo, gref, galt):
    """Compute VCF-style anchor-base alleles for an indel marker.

    VCF requires that every record shares at least one reference base between
    REF and ALT.  For indels the anchor base at position mapinfo-1 is prepended:

      deletion  (gref != '', galt == ''): pos=mapinfo-1, REF=anchor+gref, ALT=anchor
      insertion (gref == '', galt != ''): pos=mapinfo-1, REF=anchor,      ALT=anchor+galt

    fasta   : open pysam.FastaFile
    chrom   : chromosome name
    mapinfo : 1-based position of the first base of gref (or insertion site)

    Returns (vcf_pos, vcf_ref, vcf_alt).
    """
    try:
        anchor = fasta.fetch(chrom, mapinfo - 2, mapinfo - 1).upper()
        if not anchor:
            anchor = "N"
    except (ValueError, KeyError):
        anchor = "N"
    vcf_pos = mapinfo - 1
    vcf_ref = anchor + gref
    vcf_alt = anchor + galt
    return vcf_pos, vcf_ref, vcf_alt


# ── EXCLUDE INDELS FILTER ────────────────────────────────────────────────────

def apply_exclude_indels_filter(df):
    """Remove rows where _gref or _galt is empty string (indel markers).

    Called when --exclude-indels is set.  SNPs have both alleles as single
    non-empty characters; indels have one empty-string allele.
    """
    return df[(df["_gref"] != "") & (df["_galt"] != "")]


# ── COORDINATE ROLE FILTER ───────────────────────────────────────────────────

def apply_coordinate_role_filter(df, assembly, role):
    """Keep rows whose anchor_{assembly} meets the required role.

    N/A is always excluded (unmapped; already removed by Stage 1, but guarded here).
    High:     topseq_n_probe only
    Moderate: topseq_n_probe + topseq_only  (default)
    Low:      topseq_n_probe + topseq_only + probe_only
    """
    col = f"anchor_{assembly}"
    allowed = {"topseq_n_probe"}
    if role in ("Moderate", "Low"):
        allowed.add("topseq_only")
    if role == "Low":
        allowed.add("probe_only")
    return df[df[col].isin(allowed)]


# ── TIE LABEL FILTER ─────────────────────────────────────────────────────────

_TIE_RESOLVED = frozenset([
    "unique", "AS_resolved", "dAS_resolved", "NM_resolved", "CoordDelta_resolved",
])
_TIE_AVOID_SCAFFOLDS = _TIE_RESOLVED | frozenset(["scaffold_resolved"])


def apply_tie_label_filter(df, assembly, label):
    """Keep rows whose tie_{assembly} meets the required label.

    tie=ambiguous is always excluded.
    unique:          unique only
    resolved:        unique + AS/dAS/NM/CoordDelta_resolved  (default)
    avoid_scaffolds: resolved + scaffold_resolved
    """
    col = f"tie_{assembly}"
    if label == "unique":
        allowed = {"unique"}
    elif label == "resolved":
        allowed = _TIE_RESOLVED
    elif label == "avoid_scaffolds":
        allowed = _TIE_AVOID_SCAFFOLDS
    else:
        raise ValueError(f"Unknown tie label: {label!r}.")
    return df[df[col].isin(allowed)]


# ── REFALT CONFIDENCE FILTER ──────────────────────────────────────────────────

_REFALT_HIGH     = frozenset(["NM_match", "NM_validated"])
_REFALT_MODERATE = _REFALT_HIGH | frozenset(["NM_N/A", "NM_tied"])
_REFALT_LOW      = _REFALT_MODERATE | frozenset(["NM_only", "NM_unmatch", "NM_corrected"])


def apply_refalt_conf_filter(df, assembly, conf):
    """Keep rows whose RefAltMethodAgreement_{assembly} meets the required confidence.

    NM_mismatch and ambiguous are always excluded.
    High:     NM_match, NM_validated
    Moderate: High + NM_N/A, NM_tied  (default)
    Low:      Moderate + NM_only, NM_unmatch, NM_corrected
    """
    col = f"RefAltMethodAgreement_{assembly}"
    if conf == "High":
        allowed = _REFALT_HIGH
    elif conf == "Moderate":
        allowed = _REFALT_MODERATE
    elif conf == "Low":
        allowed = _REFALT_LOW
    else:
        raise ValueError(f"Unknown refalt conf: {conf!r}.")
    return df[df[col].isin(allowed)]


# ── 3D TABLE FORMATTER ───────────────────────────────────────────────────────

def format_three_d_table(three_d):
    """Format a 3-Dimension Summary (anchor × tie × RefAlt bucket) as a string.

    three_d: dict mapping (anchor, tie) → {"NM_*": int, "ambiguous": int, "N/A": int}
    """
    ANCHOR_ORDER = ["topseq_n_probe", "topseq_only", "probe_only", "N/A"]
    TIE_ORDER    = ["unique", "AS_resolved", "dAS_resolved", "NM_resolved",
                    "CoordDelta_resolved", "scaffold_resolved", "ambiguous", "N/A"]
    W = 70

    lines = [
        "═" * W,
        "3-Dimension Summary  (anchor × tie × Ref/Alt outcome)  — final markers",
        f"  {'anchor / tie':<28} {'NM_*(Chr≠0)':>10} {'amb(Chr=0)':>10}"
        f" {'N/A(Chr=0)':>10} {'Total':>8}",
        "  " + "─" * 68,
    ]

    grand = {"NM_*": 0, "ambiguous": 0, "N/A": 0}
    for anchor in ANCHOR_ORDER:
        anchor_data = {t: d for (a, t), d in three_d.items() if a == anchor}
        if sum(v for d in anchor_data.values() for v in d.values()) == 0:
            continue
        lines.append(f"  anchor={anchor}")
        for tie in TIE_ORDER:
            d = anchor_data.get(tie, {})
            nm, amb, na = d.get("NM_*", 0), d.get("ambiguous", 0), d.get("N/A", 0)
            if nm + amb + na == 0:
                continue
            lines.append(
                f"    tie={tie:<24} {nm:>10,} {amb:>10,} {na:>10,} {nm+amb+na:>8,}"
            )
            grand["NM_*"] += nm
            grand["ambiguous"] += amb
            grand["N/A"] += na

    total = sum(grand.values())
    lines += [
        "  " + "─" * 68,
        f"  {'Total':<28} {grand['NM_*']:>10,} {grand['ambiguous']:>10,}"
        f" {grand['N/A']:>10,} {total:>8,}",
    ]
    return "\n".join(lines)


# ── MAIN ─────────────────────────────────────────────────────────────────────

def run_qc(args):
    assembly   = args.assembly
    col_chr    = f"Chr_{assembly}"
    col_pos    = f"MapInfo_{assembly}"
    col_strand = f"Strand_{assembly}"
    col_ref    = f"Ref_{assembly}"
    col_alt    = f"Alt_{assembly}"

    out_dir = os.path.abspath(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)
    temp_dir = os.path.abspath(args.temp_dir) if args.temp_dir else out_dir
    os.makedirs(temp_dir, exist_ok=True)
    assessment_dir = os.path.join(out_dir, "remap_assessment")
    os.makedirs(assessment_dir, exist_ok=True)

    # Derive prefix from input filename if not given
    prefix = args.prefix or os.path.splitext(os.path.basename(args.input))[0]
    prefix = re.sub(r"_remapped$", "", prefix)

    qc_stats = {}

    # ── Load remapped manifest ───────────────────────────────────────────────
    print(f"[qc] Loading remapped manifest: {args.input}")
    df = pd.read_csv(args.input, low_memory=False, dtype={col_chr: str, "Chr": str})
    qc_stats["Input markers"] = len(df)
    print(f"[qc] {len(df):,} markers loaded.")

    # ── Benchmark against known-assembly markers ─────────────────────────────
    benchmark_known_assembly(df, assembly, assessment_dir)

    # ── MAPQ histograms (all markers) ────────────────────────────────────────
    write_mapq_histo(df["MAPQ_TopGenomicSeq"].dropna(), 2,
                     os.path.join(assessment_dir, "MAPQ_TopGenomicSeq.histo"))
    write_mapq_histo(df["MAPQ_Probe"].dropna(), 2,
                     os.path.join(assessment_dir, "MAPQ_Probe.histo"))

    # ── Filter 1: Strand=N/A (unmapped + ambiguous) ─────────────────────────
    df_mapped = df[df[col_strand].isin(["+", "-"])].copy()
    qc_stats["After Strand=N/A filter (unmapped + ambiguous)"] = len(df_mapped)
    print(f"[qc] After Strand=N/A filter (unmapped + ambiguous): {len(df_mapped):,} ({len(df) - len(df_mapped):,} removed)")

    # ── Filter 2: MAPQ filter ────────────────────────────────────────────────
    # NaN MAPQ_TopGenomicSeq means probe_only (no TopSeq alignment) — exempt from this filter.
    ts_mapq = df_mapped["MAPQ_TopGenomicSeq"]
    if args.mapq_topseq > 0:
        ts_fail = ts_mapq.notna() & (ts_mapq < args.mapq_topseq)
        df_topseq_filtered = df_mapped[~ts_fail]
    else:
        df_topseq_filtered = df_mapped
    df_mapq = apply_probe_mapq_filter(df_topseq_filtered, args.mapq_probe).copy()
    qc_stats["After MAPQ filter"] = len(df_mapq)
    print(f"[qc] After MAPQ filter (TopSeq>={args.mapq_topseq}): {len(df_mapq):,} ({len(df_mapped) - len(df_mapq):,} removed)")

    # ── Filter 2.5: CoordDelta filter ────────────────────────────────────────
    col_coord_delta  = f"CoordDelta_{assembly}"
    col_anchor       = f"anchor_{assembly}"
    if args.coord_delta >= 0:
        if col_coord_delta not in df_mapq.columns or col_anchor not in df_mapq.columns:
            print(f"[qc] WARNING: --coord-delta requested but {col_coord_delta!r} or "
                  f"{col_anchor!r} column not found in input. Skipping filter. "
                  "Re-run remap_manifest.py to generate these columns.")
            df_coord = df_mapq
        else:
            # Remove markers where CoordDelta > threshold.
            # CoordDelta = -1 means CIGAR coord was unavailable (SNP in soft-clipped region);
            # these are probe-derived coordinates with no cross-validation — they pass through
            # (−1 is not > any threshold ≥ 0).
            # topseq_only markers also carry CoordDelta = -1 but have no probe at all;
            # they are explicitly removed whenever the filter is active.
            exceeds_delta  = df_mapq[col_coord_delta] > args.coord_delta
            is_topseq_only = df_mapq[col_anchor] == "topseq_only"
            df_coord = df_mapq[~exceeds_delta & ~is_topseq_only].copy()
            n_removed = len(df_mapq) - len(df_coord)
            n_delta   = exceeds_delta.sum()
            n_ts_only = (is_topseq_only & ~exceeds_delta).sum()
            print(f"[qc] After CoordDelta filter (delta<={args.coord_delta}): "
                  f"{len(df_coord):,} ({n_removed:,} removed: "
                  f"{n_delta:,} CoordDelta>{args.coord_delta}, {n_ts_only:,} topseq_only)")
    else:
        df_coord = df_mapq
    if args.coord_delta >= 0:
        qc_stats[f"After CoordDelta filter (delta<={args.coord_delta})"] = len(df_coord)

    # ── Filter 2.6: Strand agreement filter ─────────────────────────────────
    col_strand_agree = f"StrandAgreementAsExpected_{assembly}"
    if args.require_strand_agreement:
        if col_strand_agree not in df_coord.columns:
            print(f"[qc] WARNING: --require-strand-agreement requested but {col_strand_agree!r} "
                  "column not found in input. Skipping filter. "
                  "Re-run remap_manifest.py to generate this column.")
        else:
            df_before = df_coord
            df_coord = apply_strand_agreement_filter(df_coord, assembly).copy()
            n_removed = len(df_before) - len(df_coord)
            print(f"[qc] After strand agreement filter: {len(df_coord):,} ({n_removed:,} removed)")
            qc_stats["After strand agreement filter"] = len(df_coord)

    # ── Generate VCF of all remapped positions ───────────────────────────────
    print("[qc] Generating VCF position template...")
    pos_vcf  = os.path.join(temp_dir, "_pos.vcf")
    ref_vcf  = os.path.join(temp_dir, "_ref.vcf")
    build_pos_vcf(df_coord, args.vcf_contigs, col_chr, col_pos, pos_vcf)

    print("[qc] Extracting reference alleles with bcftools...")
    ref_alleles = extract_ref_alleles(pos_vcf, args.reference, ref_vcf)

    # ── Strand-normalise remapped alleles (to + strand) ─────────────────────
    df_coord = df_coord.copy()
    df_coord["_gref"] = df_coord.apply(
        lambda r: strand_normalize(str(r[col_ref]) if pd.notna(r[col_ref]) else "", r[col_strand]), axis=1)
    df_coord["_galt"] = df_coord.apply(
        lambda r: strand_normalize(str(r[col_alt]) if pd.notna(r[col_alt]) else "", r[col_strand]), axis=1)
    df_coord["_genome_ref"] = df_coord["Name"].map(ref_alleles)

    # ── Auto-correct swapped Ref/Alt assignments ────────────────────────────
    swap_mask = (
        (df_coord["_gref"] != df_coord["_genome_ref"]) &
        (df_coord["_galt"] == df_coord["_genome_ref"])
    )
    if swap_mask.any():
        n_swapped = swap_mask.sum()
        print(f"[qc] Auto-correcting {n_swapped:,} swapped Ref/Alt assignments (Alt matched genome Ref).")
        df_coord.loc[swap_mask, ["_gref", "_galt"]] = (
            df_coord.loc[swap_mask, ["_galt", "_gref"]].values
        )
        df_coord.loc[swap_mask, [col_ref, col_alt]] = (
            df_coord.loc[swap_mask, [col_alt, col_ref]].values
        )

    # ── Filter 3: Design conflict (remapped Ref must match genome Ref) ───────
    # For SNPs: _gref (single base, + strand) must equal _genome_ref from bcftools.
    # For indels with deletion ref (_gref != ''): use pysam to fetch the full
    #   deletion sequence from the reference and compare to _gref.
    # For insertions (_gref == ''): no reference sequence to verify; pass through.
    ref_fasta = pysam.FastaFile(args.reference)
    try:
        snp_mask   = (df_coord["_gref"] != "") & (df_coord["_galt"] != "")
        indel_mask = (df_coord["_gref"] == "") | (df_coord["_galt"] == "")

        # SNP design-conflict check (original logic)
        snp_pass = snp_mask & (df_coord["_gref"] == df_coord["_genome_ref"])

        # Indel design-conflict check: use RefAltMethodAgreement computed during remapping.
        # 'NM_mismatch' means the NM-determined Ref allele did not match the genome sequence.
        # If the column is absent (old CSV), fall back to the original pysam check.
        col_refalt_agree = f"RefAltMethodAgreement_{assembly}"
        if col_refalt_agree in df_coord.columns:
            indel_pass = indel_mask & (df_coord[col_refalt_agree] != "NM_mismatch")
        else:
            def _indel_passes(row):
                return check_deletion_ref_match(
                    ref_fasta, row[col_chr], int(row[col_pos]), row["_gref"]
                )
            indel_pass = indel_mask & df_coord.apply(_indel_passes, axis=1)

        matching = df_coord[snp_pass | indel_pass].copy()
    finally:
        ref_fasta.close()

    qc_stats["After design conflict filter"] = len(matching)
    print(f"[qc] After design conflict filter: {len(matching):,} ({len(df_coord) - len(matching):,} removed)")

    # ── Filter 3.5 (optional): Exclude indel markers ─────────────────────────
    if args.exclude_indels:
        before_excl = len(matching)
        matching = apply_exclude_indels_filter(matching).copy()
        n_excl = before_excl - len(matching)
        print(f"[qc] After exclude-indels filter: {len(matching):,} ({n_excl:,} indels removed)")
        qc_stats["After exclude-indels filter"] = len(matching)

    # Write _matchingSNPs.vcf
    match_vcf = os.path.join(out_dir, "_matchingSNPs.vcf")
    ref_fasta2 = pysam.FastaFile(args.reference)
    try:
        with open(match_vcf, "w") as f:
            with open(args.vcf_contigs) as vc:
                f.write("##fileformat=VCFv4.3\n")
                f.write(vc.read())
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for _, row in matching.iterrows():
                gref = row["_gref"]
                galt = row["_galt"]
                pos  = int(row[col_pos])
                if gref == "" or galt == "":
                    # Indel: use anchor-base VCF encoding
                    vcf_pos, vcf_ref, vcf_alt = make_anchor_alleles(
                        ref_fasta2, row[col_chr], pos, gref, galt
                    )
                else:
                    vcf_pos, vcf_ref, vcf_alt = pos, gref, galt
                f.write(f"{row[col_chr]}\t{vcf_pos}\t{row['Name']}\t{vcf_ref}\t{vcf_alt}\t.\t.\t.\n")
    finally:
        ref_fasta2.close()

    # ── Filter 4: Polymorphic positions ─────────────────────────────────────
    # A position is polymorphic if multiple markers at the same Chr:Pos have different Ref/Alt
    pos_allele_counts = (
        matching.assign(_allele_pair=matching["_gref"] + "," + matching["_galt"])
        .groupby([col_chr, col_pos])["_allele_pair"]
        .nunique()
    )
    polymorphic = pos_allele_counts[pos_allele_counts > 1].reset_index()
    poly_set = set(zip(polymorphic[col_chr], polymorphic[col_pos]))

    poly_path = os.path.join(out_dir, "polymorphic_positions.txt")
    with open(poly_path, "w") as f:
        for chr_, pos in poly_set:
            f.write(f"{chr_},{pos}\n")

    df_binary = matching[~matching.apply(
        lambda r: (r[col_chr], r[col_pos]) in poly_set, axis=1
    )].copy()
    qc_stats["After polymorphic filter"] = len(df_binary)
    print(f"[qc] After polymorphic filter: {len(df_binary):,} ({len(matching) - len(df_binary):,} removed)")

    binary_vcf = os.path.join(out_dir, "_matchingSNPs_binary.vcf")
    with open(binary_vcf, "w") as f:
        with open(match_vcf) as src:
            for line in src:
                if line.startswith("#"):
                    f.write(line)
                    continue
                name = line.split("\t")[2]
                if name in df_binary["Name"].values:
                    f.write(line)

    # ── Filter 5: Consistency check ─────────────────────────────────────────
    topseq_sam = args.topseq_sam or os.path.join(temp_dir, "temp_topseq.sam")
    probe_sam  = args.probe_sam  or os.path.join(temp_dir, "temp_probe.sam")

    if os.path.exists(topseq_sam) and os.path.exists(probe_sam):
        print("[qc] Running consistency check against SAM files...")
        inconsistent = get_consistent_snps(topseq_sam, probe_sam)

        # MAPQ histograms for inconsistent markers
        df_incons = df[df["Name"].isin(inconsistent)]
        write_mapq_histo(df_incons["MAPQ_TopGenomicSeq"].dropna(), 2,
                         os.path.join(assessment_dir, "MAPQ_TopGenomicSeq_inconsistent_remapped.histo"))
        write_mapq_histo(df_incons["MAPQ_Probe"].dropna(), 2,
                         os.path.join(assessment_dir, "MAPQ_Probe_inconsistent_remapped.histo"))
        df_incons.to_csv(os.path.join(out_dir, "inconsistent_remapped.csv"), index=False)

        df_consistent = df_binary[~df_binary["Name"].isin(inconsistent)].copy()
    else:
        print(f"[qc] WARNING: SAM files not found ({topseq_sam}, {probe_sam}). "
              "Skipping consistency filter. Pass --topseq-sam and --probe-sam to enable.")
        inconsistent = set()
        df_consistent = df_binary.copy()

    qc_stats["After consistency filter"] = len(df_consistent)
    print(f"[qc] After consistency filter: {len(df_consistent):,} ({len(df_binary) - len(df_consistent):,} removed)")

    # Write _matchingSNPs_binary_consistantMapping.vcf
    consist_vcf = os.path.join(out_dir, "_matchingSNPs_binary_consistantMapping.vcf")
    consist_names = set(df_consistent["Name"])
    with open(binary_vcf) as src, open(consist_vcf, "w") as f:
        for line in src:
            if line.startswith("#"):
                f.write(line)
                continue
            if line.split("\t")[2] in consist_names:
                f.write(line)

    # ── PLINK BIM file ───────────────────────────────────────────────────────
    # For indel markers, apply anchor-base encoding so PLINK sees VCF-style alleles.
    bim_path = os.path.join(out_dir, f"{prefix}_remapped_{assembly}.bim")
    bim = df_consistent[[col_chr, "Name", col_pos, "_gref", "_galt"]].copy()
    ref_fasta3 = pysam.FastaFile(args.reference)
    try:
        def _bim_row_alleles(row):
            gref = row["_gref"]
            galt = row["_galt"]
            if gref == "" or galt == "":
                _, vcf_ref, vcf_alt = make_anchor_alleles(
                    ref_fasta3, row[col_chr], int(row[col_pos]), gref, galt
                )
                return vcf_ref, vcf_alt
            return gref, galt
        bim_alleles = bim.apply(_bim_row_alleles, axis=1, result_type="expand")
        bim["_gref"] = bim_alleles[0]
        bim["_galt"] = bim_alleles[1]
    finally:
        ref_fasta3.close()

    bim.insert(2, "cM", 0)
    bim = bim.sort_values([col_chr, col_pos])
    bim.to_csv(bim_path, sep="\t", header=False, index=False)
    print(f"[qc] BIM written: {bim_path}")

    # ── Ambiguous SNP count ──────────────────────────────────────────────────
    ambig_pairs = {frozenset(["A","T"]), frozenset(["C","G"])}
    ambiguous = bim[bim.apply(
        lambda r: frozenset([str(r["_gref"]), str(r["_galt"])]) in ambig_pairs, axis=1
    )]
    qc_stats["Ambiguous SNPs (A/T or C/G)"] = len(ambiguous)

    # ── Final map file ───────────────────────────────────────────────────────
    map_path = os.path.join(out_dir, f"matchingSNPs_binary_consistantMapping.{assembly}_map")
    print("[qc] Building final map file...")
    errors = build_final_map(df_consistent, assembly, map_path)
    if errors:
        print(f"[qc] WARNING: {errors} markers could not be assigned to map file (written as 'Error' lines).")
    else:
        print(f"[qc] Final map file written with 0 errors: {map_path}")

    qc_stats["Final markers"] = len(df_consistent)

    # ── QC Report ────────────────────────────────────────────────────────────
    report_path = os.path.join(out_dir, "QC_Report.txt")
    with open(report_path, "w") as f:
        f.write(f"QC Report — assembly: {assembly}\n")
        f.write(f"Input: {args.input}\n")
        f.write("-" * 60 + "\n")
        prev = None
        for stage, count in qc_stats.items():
            if prev is not None:
                removed = prev - count if isinstance(count, int) and isinstance(prev, int) else ""
                removed_str = f"  (-{prev - count:,})" if removed != "" else ""
            else:
                removed_str = ""
            f.write(f"{stage:<40} {count:>8,}{removed_str}\n")
            if isinstance(count, int):
                prev = count
    print(f"[qc] QC report written: {report_path}")

    # Print summary
    print("\n--- QC Summary ---")
    with open(report_path) as f:
        print(f.read())

    print("[qc] Done.")


if __name__ == "__main__":
    run_qc(parse_args())
