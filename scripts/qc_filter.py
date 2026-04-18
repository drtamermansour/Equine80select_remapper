"""
qc_filter.py — QC filtering, allele decision, VCF/BIM/map generation.

Consumes the remapped manifest produced by ``remap_manifest.py`` and applies an
11-stage filter cascade, then writes downstream output files. Per-stage counts
are recorded in ``QC_Report.txt`` and a per-marker trace CSV annotates each
input marker with the first filter stage (if any) that rejected it
(``WhyFiltered_{assembly}`` column). See ``WHY_FILTERED_LABELS`` below for the
canonical stage order.

Outputs under ``--output-dir``:

  {prefix}_allele_map_{assembly}.tsv                       Manifest <-> genome allele crosswalk (main output)
  {prefix}_remapped_{assembly}.bim                         PLINK BIM
  {prefix}_remapped_{assembly}.vcf                         Final filtered VCF
  {prefix}_remapped_{assembly}_traced.csv                  Full input + WhyFiltered column
  QC_Report.txt                                            Per-stage counts
  diagnostics/                                             MAPQ histograms + known-assembly benchmark

Usage:
  python scripts/qc_filter.py \\
      -i remapped.csv \\
      -r reference.fa \\
      -v vcf_contigs.txt \\
      -a equCab3 \\
      -o output_dir/ \\
      [--coordinate-role Moderate] [--tie-label resolved] [--refalt-conf Moderate] \\
      [--mapq-topseq 30] [--mapq-probe 0] [--coord-delta -1] \\
      [--keep-indels] [--keep-polymorphic] [--keep-ambiguous] \\
      [--temp-dir /tmp/remap] [--prefix Equine80select]
"""

import argparse
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd
import pysam

from _strand_utils import strand_normalize, complement


# ── Per-marker filter tracing (WhyFiltered_{assembly} column) ────────────────

WHY_FILTERED_LABELS = [
    "stage_1_failed_markers",
    "stage_2_design_conflict",
    "stage_3_coordinate_role",
    "stage_4_tie_label",
    "stage_5_refalt_conf",
    "stage_6_mapq_topseq",
    "stage_7_mapq_probe",
    "stage_8_coord_delta",
    "stage_9_indel_excluded",
    "stage_10_polymorphic",
    "stage_11_ambiguous_snp",
]


def _tag_removed(why: pd.Series, before_idx, after_idx, label: str) -> None:
    """Mark markers that survived into *before_idx* but not *after_idx* with *label*.

    First-rejection-wins: markers already tagged by an earlier stage keep that tag.
    Mutates *why* in place; caller owns the series.
    """
    removed = pd.Index(before_idx).difference(pd.Index(after_idx))
    if len(removed) == 0:
        return
    untagged = removed[why.loc[removed] == ""]
    why.loc[untagged] = label


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
    p.add_argument(
        "--keep-ambiguous", action="store_true", default=False,
        help="Keep ambiguous (A/T or C/G) SNPs in outputs (default: False = removed)."
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


# ── FINAL MAP FILE ───────────────────────────────────────────────────────────

def build_final_map(df_final, assembly, map_path):
    """
    Writes {prefix}_allele_map_{assembly}.tsv with a header row and columns:
      chr  pos  snp_id  manifest_alleles  genomic_alleles  manifest_ref  genomic_ref  decision

    manifest_alleles are from the manifest SNP column (e.g. A,G).
    genomic_alleles are the + strand remapped alleles.
    decision ('as_is' or 'complement') is inferred by matching manifest alleles to genomic alleles —
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
        f.write("chr\tpos\tsnp_id\tmanifest_alleles\tgenomic_alleles\tmanifest_ref\tgenomic_ref\tdecision\n")
        if lines:
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
    Writes a mismatch report to diagnostics/{assembly}.mismatches.
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


# ── EXCLUDE AMBIGUOUS-SNPS FILTER ─────────────────────────────────────────────

_AMBIG_PAIRS = frozenset({frozenset({"A", "T"}), frozenset({"C", "G"})})


def polymorphic_positions(df, col_chr, col_pos):
    """Return the set of (chr, pos) tuples where multiple markers emit conflicting
    Ref/Alt assignments at the same coordinate.

    Used by Stage 10 both to remove polymorphic markers and (when `--keep-polymorphic`
    is set) to compute the hypothetical count of markers that *would* have been
    removed if the filter had been applied.
    """
    if len(df) == 0:
        return set()
    counts = (
        df.assign(_allele_pair=df["_gref"] + "," + df["_galt"])
          .groupby([col_chr, col_pos])["_allele_pair"]
          .nunique()
    )
    polymorphic = counts[counts > 1].reset_index()
    return set(zip(polymorphic[col_chr], polymorphic[col_pos]))


def apply_exclude_ambiguous_snps_filter(df):
    """Remove SNPs whose {_gref, _galt} pair is ambiguous after strand-normalisation.

    Ambiguous pairs ({A,T} and {C,G}) equal their own reverse-complement, so
    the alleles alone cannot resolve which strand the variant lives on —
    downstream imputation and GWAS typically exclude these.

    Only SNPs trigger: indels have an empty-string allele and never match
    a two-base ambiguous set.
    """
    def _is_ambig(row):
        return frozenset({str(row["_gref"]), str(row["_galt"])}) in _AMBIG_PAIRS
    if len(df) == 0:
        return df
    return df[~df.apply(_is_ambig, axis=1)]


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

    tie=locus_unresolved is always excluded.
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

    NM_mismatch and refalt_unresolved are always excluded.
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

    three_d: dict mapping (anchor, tie) → {"NM_*": int, "refalt_unresolved": int, "N/A": int}
    """
    ANCHOR_ORDER = ["topseq_n_probe", "topseq_only", "probe_only", "N/A"]
    TIE_ORDER    = ["unique", "AS_resolved", "dAS_resolved", "NM_resolved",
                    "CoordDelta_resolved", "scaffold_resolved", "locus_unresolved", "N/A"]
    W = 70

    lines = [
        "═" * W,
        "3-Dimension Summary  (anchor × tie × Ref/Alt outcome)  — final markers",
        f"  {'anchor / tie':<28} {'NM_*(Chr≠0)':>10} {'unresolved(Chr=0)':>17}"
        f" {'not_attempted(Chr=0)':>20} {'Total':>8}",
        "  " + "─" * 86,
    ]

    grand = {"NM_*": 0, "refalt_unresolved": 0, "N/A": 0}
    for anchor in ANCHOR_ORDER:
        anchor_data = {t: d for (a, t), d in three_d.items() if a == anchor}
        if sum(v for d in anchor_data.values() for v in d.values()) == 0:
            continue
        lines.append(f"  anchor={anchor}")
        for tie in TIE_ORDER:
            d = anchor_data.get(tie, {})
            nm, amb, na = d.get("NM_*", 0), d.get("refalt_unresolved", 0), d.get("N/A", 0)
            if nm + amb + na == 0:
                continue
            lines.append(
                f"    tie={tie:<24} {nm:>10,} {amb:>17,} {na:>20,} {nm+amb+na:>8,}"
            )
            grand["NM_*"] += nm
            grand["refalt_unresolved"] += amb
            grand["N/A"] += na

    total = sum(grand.values())
    lines += [
        "  " + "─" * 86,
        f"  {'Total':<28} {grand['NM_*']:>10,} {grand['refalt_unresolved']:>17,}"
        f" {grand['N/A']:>20,} {total:>8,}",
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
    assessment_dir = os.path.join(out_dir, "diagnostics")
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

    # Per-marker filter trace (populated at each stage; empty for passing markers)
    col_why = f"WhyFiltered_{assembly}"
    why_filtered = pd.Series("", index=df.index, dtype="object")

    # ── Benchmark against known-assembly markers ─────────────────────────────
    benchmark_known_assembly(df, assembly, assessment_dir)

    # ── MAPQ histograms (all markers) ────────────────────────────────────────
    write_mapq_histo(df["MAPQ_TopGenomicSeq"].dropna(), 2,
                     os.path.join(assessment_dir, "MAPQ_TopGenomicSeq.histo"))
    write_mapq_histo(df["MAPQ_Probe"].dropna(), 2,
                     os.path.join(assessment_dir, "MAPQ_Probe.histo"))

    # ── Stage 1: Failed markers (Strand=N/A — unmapped + locus_unresolved) ──
    df_mapped = df[df[col_strand].isin(["+", "-"])].copy()
    _tag_removed(why_filtered, df.index, df_mapped.index, "stage_1_failed_markers")
    qc_stats["Failed markers (unmapped + locus_unresolved)"] = len(df_mapped)
    print(f"[qc] Stage 1 — Failed markers removed: {len(df) - len(df_mapped):,}; remaining: {len(df_mapped):,}")

    # ── VCF generation + strand normalisation (needed by Stage 2) ───────────
    print("[qc] Generating VCF position template...")
    pos_vcf  = os.path.join(temp_dir, "_pos.vcf")
    ref_vcf  = os.path.join(temp_dir, "_ref.vcf")
    build_pos_vcf(df_mapped, args.vcf_contigs, col_chr, col_pos, pos_vcf)

    print("[qc] Extracting reference alleles with bcftools...")
    ref_alleles = extract_ref_alleles(pos_vcf, args.reference, ref_vcf)

    df_mapped["_gref"] = df_mapped.apply(
        lambda r: strand_normalize(str(r[col_ref]) if pd.notna(r[col_ref]) else "", r[col_strand]), axis=1)
    df_mapped["_galt"] = df_mapped.apply(
        lambda r: strand_normalize(str(r[col_alt]) if pd.notna(r[col_alt]) else "", r[col_strand]), axis=1)
    df_mapped["_genome_ref"] = df_mapped["Name"].map(ref_alleles)

    # ── Auto-correct swapped Ref/Alt assignments ────────────────────────────
    swap_mask = (
        (df_mapped["_gref"] != df_mapped["_genome_ref"]) &
        (df_mapped["_galt"] == df_mapped["_genome_ref"])
    )
    if swap_mask.any():
        n_swapped = swap_mask.sum()
        print(f"[qc] Auto-correcting {n_swapped:,} swapped Ref/Alt assignments (Alt matched genome Ref).")
        df_mapped.loc[swap_mask, ["_gref", "_galt"]] = (
            df_mapped.loc[swap_mask, ["_galt", "_gref"]].values
        )
        df_mapped.loc[swap_mask, [col_ref, col_alt]] = (
            df_mapped.loc[swap_mask, [col_alt, col_ref]].values
        )

    # ── Stage 2: Design conflict (remapped Ref must match genome Ref) ────────
    # SNPs:       _gref (single base, + strand) must equal _genome_ref (from bcftools).
    # Deletions:  RefAltMethodAgreement_{assembly} != 'NM_mismatch' (already validated in remap_manifest).
    # Insertions: _gref == '' — no reference sequence to verify; pass through.
    snp_mask   = (df_mapped["_gref"] != "") & (df_mapped["_galt"] != "")
    indel_mask = (df_mapped["_gref"] == "") | (df_mapped["_galt"] == "")

    snp_pass = snp_mask & (df_mapped["_gref"] == df_mapped["_genome_ref"])

    col_refalt_agree = f"RefAltMethodAgreement_{assembly}"
    if col_refalt_agree not in df_mapped.columns:
        raise ValueError(
            f"Column {col_refalt_agree!r} is required for Stage 2. "
            f"Input CSV must be produced by the current remap_manifest.py."
        )
    indel_pass = indel_mask & (df_mapped[col_refalt_agree] != "NM_mismatch")

    df_noconflict = df_mapped[snp_pass | indel_pass].copy()

    _tag_removed(why_filtered, df_mapped.index, df_noconflict.index, "stage_2_design_conflict")
    qc_stats["After design conflict filter"] = len(df_noconflict)
    print(f"[qc] Stage 2 — Design conflict removed: {len(df_mapped) - len(df_noconflict):,}; remaining: {len(df_noconflict):,}")

    # ── Stage 3: Coordinate role ─────────────────────────────────────────────
    col_anchor = f"anchor_{assembly}"
    if col_anchor in df_noconflict.columns:
        df_coord_role = apply_coordinate_role_filter(df_noconflict, assembly, args.coordinate_role).copy()
        n_removed = len(df_noconflict) - len(df_coord_role)
        print(f"[qc] Stage 3 — Coordinate role ({args.coordinate_role}): {n_removed:,} removed; {len(df_coord_role):,} remaining")
    else:
        print(f"[qc] WARNING: {col_anchor!r} column not found. Skipping coordinate role filter.")
        df_coord_role = df_noconflict
    _tag_removed(why_filtered, df_noconflict.index, df_coord_role.index, "stage_3_coordinate_role")
    qc_stats[f"After coordinate role ({args.coordinate_role})"] = len(df_coord_role)

    # ── Stage 4: Tie label ────────────────────────────────────────────────────
    col_tie = f"tie_{assembly}"
    if col_tie in df_coord_role.columns:
        df_tie = apply_tie_label_filter(df_coord_role, assembly, args.tie_label).copy()
        n_removed = len(df_coord_role) - len(df_tie)
        print(f"[qc] Stage 4 — Tie label ({args.tie_label}): {n_removed:,} removed; {len(df_tie):,} remaining")
    else:
        print(f"[qc] WARNING: {col_tie!r} column not found. Skipping tie label filter.")
        df_tie = df_coord_role
    _tag_removed(why_filtered, df_coord_role.index, df_tie.index, "stage_4_tie_label")
    qc_stats[f"After tie label ({args.tie_label})"] = len(df_tie)

    # ── Stage 5: RefAlt confidence ────────────────────────────────────────────
    if col_refalt_agree in df_tie.columns:
        df_refalt = apply_refalt_conf_filter(df_tie, assembly, args.refalt_conf).copy()
        n_removed = len(df_tie) - len(df_refalt)
        print(f"[qc] Stage 5 — RefAlt confidence ({args.refalt_conf}): {n_removed:,} removed; {len(df_refalt):,} remaining")
    else:
        print(f"[qc] WARNING: {col_refalt_agree!r} column not found. Skipping RefAlt confidence filter.")
        df_refalt = df_tie
    _tag_removed(why_filtered, df_tie.index, df_refalt.index, "stage_5_refalt_conf")
    qc_stats[f"After RefAlt confidence ({args.refalt_conf})"] = len(df_refalt)

    # ── Stage 6: MAPQ_TopGenomicSeq (probe_only exempt via NaN) ──────────────
    ts_mapq = df_refalt["MAPQ_TopGenomicSeq"]
    if args.mapq_topseq > 0:
        ts_fail = ts_mapq.notna() & (ts_mapq < args.mapq_topseq)
        df_mapq_ts = df_refalt[~ts_fail].copy()
        _tag_removed(why_filtered, df_refalt.index, df_mapq_ts.index, "stage_6_mapq_topseq")
        n_removed = len(df_refalt) - len(df_mapq_ts)
        print(f"[qc] Stage 6 — MAPQ_TopGenomicSeq (>={args.mapq_topseq}): {n_removed:,} removed; {len(df_mapq_ts):,} remaining")
        qc_stats[f"After MAPQ_TopGenomicSeq (>={args.mapq_topseq})"] = len(df_mapq_ts)
    else:
        df_mapq_ts = df_refalt
        print("[qc] Stage 6 — MAPQ_TopGenomicSeq skipped (--mapq-topseq 0).")
        qc_stats["Stage 6 skipped — MAPQ_TopGenomicSeq (--mapq-topseq 0)"] = len(df_mapq_ts)

    # ── Stage 7: MAPQ_Probe (topseq_only exempt via NaN) ─────────────────────
    if args.mapq_probe > 0:
        df_mapq = apply_probe_mapq_filter(df_mapq_ts, args.mapq_probe).copy()
        _tag_removed(why_filtered, df_mapq_ts.index, df_mapq.index, "stage_7_mapq_probe")
        n_removed = len(df_mapq_ts) - len(df_mapq)
        print(f"[qc] Stage 7 — MAPQ_Probe (>={args.mapq_probe}): {n_removed:,} removed; {len(df_mapq):,} remaining")
        qc_stats[f"After MAPQ_Probe (>={args.mapq_probe})"] = len(df_mapq)
    else:
        df_mapq = df_mapq_ts
        print(f"[qc] Stage 7 — MAPQ_Probe skipped (--mapq-probe {args.mapq_probe}).")
        qc_stats[f"Stage 7 skipped — MAPQ_Probe (--mapq-probe {args.mapq_probe})"] = len(df_mapq)

    # ── Stage 8: CoordDelta ───────────────────────────────────────────────────
    # topseq_only and probe_only have CoordDelta=-1 (no CIGAR coord) and pass through.
    # Only markers with real CoordDelta > threshold are excluded.
    col_coord_delta = f"CoordDelta_{assembly}"
    if args.coord_delta >= 0:
        if col_coord_delta in df_mapq.columns:
            exceeds = df_mapq[col_coord_delta] > args.coord_delta
            df_coord = df_mapq[~exceeds].copy()
            n_removed = len(df_mapq) - len(df_coord)
            print(f"[qc] Stage 8 — CoordDelta (delta<={args.coord_delta}): "
                  f"{n_removed:,} removed (CoordDelta>{args.coord_delta} only); {len(df_coord):,} remaining")
        else:
            print(f"[qc] WARNING: {col_coord_delta!r} column not found. Skipping CoordDelta filter.")
            df_coord = df_mapq
        qc_stats[f"After CoordDelta filter (delta<={args.coord_delta})"] = len(df_coord)
    else:
        df_coord = df_mapq
        print(f"[qc] Stage 8 — CoordDelta skipped (--coord-delta {args.coord_delta}).")
        qc_stats[f"Stage 8 skipped — CoordDelta (--coord-delta {args.coord_delta})"] = len(df_coord)
    _tag_removed(why_filtered, df_mapq.index, df_coord.index, "stage_8_coord_delta")

    # ── Stage 9: Indels (excluded by default; --keep-indels to include) ───────
    if not args.keep_indels:
        df_noindel = apply_exclude_indels_filter(df_coord).copy()
        n_removed = len(df_coord) - len(df_noindel)
        print(f"[qc] Stage 9 — Indels excluded: {n_removed:,} removed; {len(df_noindel):,} remaining")
        qc_stats["After indel filter"] = len(df_noindel)
    else:
        df_noindel = df_coord
        hypothetical = len(df_coord) - len(apply_exclude_indels_filter(df_coord))
        print(f"[qc] Stage 9 — Indels kept (--keep-indels); {hypothetical:,} would have been removed.")
        qc_stats[f"Stage 9 skipped — Indels (--keep-indels; would remove {hypothetical:,})"] = len(df_noindel)
    _tag_removed(why_filtered, df_coord.index, df_noindel.index, "stage_9_indel_excluded")

    # ── Stage 10: Polymorphic (removed by default; --keep-polymorphic to skip) ─
    if not args.keep_polymorphic:
        poly_set = polymorphic_positions(df_noindel, col_chr, col_pos)
        df_final = df_noindel[~df_noindel.apply(
            lambda r: (r[col_chr], r[col_pos]) in poly_set, axis=1
        )].copy()
        n_removed = len(df_noindel) - len(df_final)
        print(f"[qc] Stage 10 — Polymorphic removed: {n_removed:,}; {len(df_final):,} remaining")
        qc_stats["After polymorphic filter"] = len(df_final)
    else:
        df_final = df_noindel
        poly_set = polymorphic_positions(df_noindel, col_chr, col_pos)
        hypothetical = df_noindel.apply(
            lambda r: (r[col_chr], r[col_pos]) in poly_set, axis=1
        ).sum()
        print(f"[qc] Stage 10 — Polymorphic sites kept (--keep-polymorphic); {hypothetical:,} would have been removed.")
        qc_stats[f"Stage 10 skipped — Polymorphic (--keep-polymorphic; would remove {hypothetical:,})"] = len(df_final)
    _tag_removed(why_filtered, df_noindel.index, df_final.index, "stage_10_polymorphic")

    # ── Stage 11: Ambiguous SNPs (A/T, C/G) — excluded by default ────────────
    df_before_ambig = df_final
    if not args.keep_ambiguous:
        df_final = apply_exclude_ambiguous_snps_filter(df_before_ambig).copy()
        n_removed = len(df_before_ambig) - len(df_final)
        print(f"[qc] Stage 11 — Ambiguous SNPs excluded: {n_removed:,} removed; {len(df_final):,} remaining")
        qc_stats["After ambiguous-SNP filter"] = len(df_final)
    else:
        hypothetical = len(df_before_ambig) - len(apply_exclude_ambiguous_snps_filter(df_before_ambig))
        print(f"[qc] Stage 11 — Ambiguous SNPs kept (--keep-ambiguous); {hypothetical:,} would have been removed.")
        qc_stats[f"Stage 11 skipped — Ambiguous SNPs (--keep-ambiguous; would remove {hypothetical:,})"] = len(df_final)
    _tag_removed(why_filtered, df_before_ambig.index, df_final.index, "stage_11_ambiguous_snp")

    # ── Write full-input trace CSV with per-marker WhyFiltered column ────────
    df_trace = df.copy()
    df_trace[col_why] = why_filtered
    trace_path = os.path.join(out_dir, f"{prefix}_remapped_{assembly}_traced.csv")
    df_trace.to_csv(trace_path, index=False)
    print(f"[qc] Per-marker trace written: {trace_path}")

    # ── Final filtered VCF (post-all-filters) ────────────────────────────────
    vcf_path = os.path.join(out_dir, f"{prefix}_remapped_{assembly}.vcf")
    ref_fasta_vcf = pysam.FastaFile(args.reference)
    try:
        with open(vcf_path, "w") as f:
            with open(args.vcf_contigs) as vc:
                f.write("##fileformat=VCFv4.3\n")
                f.write(vc.read())
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for _, row in df_final.iterrows():
                gref = row["_gref"]
                galt = row["_galt"]
                pos  = int(row[col_pos])
                if gref == "" or galt == "":
                    vcf_pos, vcf_ref, vcf_alt = make_anchor_alleles(
                        ref_fasta_vcf, row[col_chr], pos, gref, galt
                    )
                else:
                    vcf_pos, vcf_ref, vcf_alt = pos, gref, galt
                f.write(f"{row[col_chr]}\t{vcf_pos}\t{row['Name']}\t{vcf_ref}\t{vcf_alt}\t.\t.\t.\n")
    finally:
        ref_fasta_vcf.close()
    print(f"[qc] VCF written: {vcf_path}")

    # ── PLINK BIM file ───────────────────────────────────────────────────────
    # For indel markers, apply anchor-base encoding so PLINK sees VCF-style alleles.
    bim_path = os.path.join(out_dir, f"{prefix}_remapped_{assembly}.bim")
    bim = df_final[[col_chr, "Name", col_pos, "_gref", "_galt"]].copy()
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

    # ── Final allele-map file (main output) ──────────────────────────────────
    map_path = os.path.join(out_dir, f"{prefix}_allele_map_{assembly}.tsv")
    print("[qc] Building allele-map file...")
    errors = build_final_map(df_final, assembly, map_path)
    if errors:
        print(f"[qc] WARNING: {errors} markers could not be assigned to allele map (written as 'Error' lines).")
    else:
        print(f"[qc] Allele-map file written with 0 errors: {map_path}")

    # ── QC Report ────────────────────────────────────────────────────────────
    # Stage 11 already surfaces the ambiguous-SNP count (either as "removed" in
    # the applied branch or as "would remove N" in the skipped branch), so a
    # separate "Ambiguous SNPs (A/T or C/G)" diagnostic line would only
    # duplicate it. Final markers is printed as a standalone row with no diff
    # (a diff against the previous stage-row count would be meaningless).
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
            f.write(f"{stage:<75} {count:>8,}{removed_str}\n")
            if isinstance(count, int):
                prev = count
        f.write(f"{'Final markers':<75} {len(df_final):>8,}\n")

    # ── 3D summary appended to QC_Report.txt ─────────────────────────────────
    _req = [f"anchor_{assembly}", f"tie_{assembly}", f"RefAltMethodAgreement_{assembly}"]
    if all(c in df_final.columns for c in _req):
        col_a_f = f"anchor_{assembly}"
        col_t_f = f"tie_{assembly}"
        col_r_f = f"RefAltMethodAgreement_{assembly}"
        three_d = {}
        for (a, t, r), cnt in df_final.groupby([col_a_f, col_t_f, col_r_f]).size().items():
            bucket = "NM_*" if r.startswith("NM_") else ("refalt_unresolved" if r == "refalt_unresolved" else "N/A")
            key = (a, t)
            if key not in three_d:
                three_d[key] = {"NM_*": 0, "refalt_unresolved": 0, "N/A": 0}
            three_d[key][bucket] += cnt
        with open(report_path, "a") as f:
            f.write("\n" + format_three_d_table(three_d) + "\n")

    print(f"[qc] QC report written: {report_path}")

    # Print summary
    print("\n--- QC Summary ---")
    with open(report_path) as f:
        print(f.read())

    print("[qc] Done.")


if __name__ == "__main__":
    run_qc(parse_args())
