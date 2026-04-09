"""
qc_filter.py — QC filtering, allele decision, VCF/BIM/map generation.

Takes the remapped manifest (output of remap_manifest.py) and applies a cascade of
quality filters, then generates all downstream output files:

  Filter cascade (with marker counts at each stage written to QC_Report.txt):
    1. Unmapped filter     — remove markers where Strand_{assembly} == 'N/A'
    2. MAPQ filter         — remove markers below MAPQ thresholds
    3. Design conflict     — keep only markers where remapped Ref matches genome Ref
    4. Polymorphic sites   — remove positions with conflicting Ref/Alt assignments
    5. Consistency filter  — remove markers with inconsistent probe/topseq alignments

  Outputs:
    allele_usage_decision.txt
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

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(seq):
    return seq.translate(COMPLEMENT)


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="QC filter and output generation for remapped manifests.")
    p.add_argument("-i", "--input",    required=True, help="Remapped manifest CSV (from remap_manifest.py)")
    p.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p.add_argument("-v", "--vcf-contigs", required=True, help="VCF contig header file")
    p.add_argument("-a", "--assembly", default="new_assembly", help="Assembly name (must match remap_manifest.py -a)")
    p.add_argument("-o", "--output-dir", default=".", help="Output directory (default: current directory)")
    p.add_argument("--mapq-topseq", type=int, default=30,
                   help="Minimum MAPQ for TopGenomicSeq alignments (default: 30)")
    p.add_argument("--mapq-probe", type=int, default=0,
                   help="Minimum MAPQ for probe alignments when >0 (default: 0 = disabled)")
    p.add_argument("--temp-dir", default=None,
                   help="Directory for intermediate files (default: output-dir)")
    p.add_argument("--prefix", default=None,
                   help="Output file prefix (default: derived from input filename)")
    p.add_argument("--topseq-sam", default=None,
                   help="Path to temp_topseq.sam from remap step (for consistency check)")
    p.add_argument("--probe-sam", default=None,
                   help="Path to temp_probe.sam from remap step (for consistency check)")
    p.add_argument("--coord-delta", type=int, default=-1,
                   help="Maximum allowed CoordDelta (|probe_coord - CIGAR_coord|). "
                        "-1 = disabled (default). Any value >= 0 removes markers with "
                        "CoordDelta > threshold and all topseq_only markers.")
    return p.parse_args()


# ── ALLELE DECISION ──────────────────────────────────────────────────────────

def allele_usage_decision(row, strand_col, assembly):
    """
    Determines whether the manifest SNP alleles should be used 'as_is' or 'complement'
    when mapping to the new assembly.

    Rules (based on IlmnStrand, SourceStrand, SourceSeq vs TopGenomicSeq, and Strand):
      The XOR logic reduces to: flip if an odd number of the three conditions are True.
        cond1: IlmnStrand != SourceStrand
        cond2: SourceSeq  != TopGenomicSeq
        cond3: Strand == '-'
      flip = cond1 XOR cond2 XOR cond3
      decision = 'complement' if flip else 'as_is'

    For indel markers ([D/I] in SNP column), the prefix 'indel_' is prepended.
    Markers with N/A strand are excluded upstream and should not reach this function.
    """
    ilmn   = row.get("IlmnStrand", "")
    src    = row.get("SourceStrand", "")
    srcseq = row.get("SourceSeq", "")
    topseq = row.get("TopGenomicSeq", "")
    strand = row.get(strand_col, "")
    snp    = row.get("SNP", "")

    cond1 = ilmn != src
    cond2 = srcseq != topseq
    cond3 = strand == "-"
    flip = cond1 ^ cond2 ^ cond3

    decision = "complement" if flip else "as_is"

    # Indel markers get a separate label (kept for tracking; excluded from VCF)
    if re.search(r"\[D/I\]|\[I/D\]", snp or ""):
        decision = "indel_" + decision

    return decision


# ── CONSISTENCY CHECK ────────────────────────────────────────────────────────

def get_consistent_snps(topseq_sam_path, probe_sam_path):
    """
    Replicates the existing consistency logic:
      For each SNP, we expect exactly 3 SAM records all mapping to the same chromosome:
        - topseq_A, topseq_B (from temp_topseq.sam, with _A/_B suffix stripped)
        - probe (from temp_probe.sam)
      SNPs where any (name, chromosome) pair has a count != 3 are flagged as inconsistent.

    This mirrors the original bash logic:
      cut -f1,3 | sort | uniq -c | awk '{if($1!=3)print}'

    Returns a set of inconsistent SNP names.
    """
    counts = Counter()

    with open(topseq_sam_path) as f:
        for line in f:
            if line.startswith("@") or line.startswith("[M"):
                continue
            cols = line.split("\t")
            flag = int(cols[1])
            if flag & 4:
                continue
            # Strip _A / _B suffix; key is (name, chromosome) to detect cross-chrom splits
            name = re.sub(r"_[AB]$", "", cols[0])
            chrom = cols[2]
            counts[(name, chrom)] += 1

    with open(probe_sam_path) as f:
        for line in f:
            if line.startswith("@") or line.startswith("[M"):
                continue
            cols = line.split("\t")
            flag = int(cols[1])
            if flag & 4:
                continue
            chrom = cols[2]
            counts[(cols[0], chrom)] += 1

    # A marker is inconsistent if any (name, chrom) pair has count != 3
    return {name for (name, chrom), cnt in counts.items() if cnt != 3}


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
        return complement(allele)
    return allele


# ── FINAL MAP FILE ───────────────────────────────────────────────────────────

def build_final_map(df_final, decisions, assembly, map_path):
    """
    Writes matchingSNPs_binary_consistantMapping.{assembly}_map with columns:
      chr  pos  snpID  SNP_alleles  genomic_alleles  SNP_ref_allele  genomic_ref_allele  decision

    SNP_alleles are from the manifest SNP column (e.g. A,G).
    genomic_alleles are the + strand remapped alleles.
    The allele_usage_decision determines which orientation the SNP alleles correspond to.
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
        decision = decisions.get(name, "")

        # Parse SNP alleles from manifest [A/G] format
        m = re.search(r"\[(.+?)/(.+?)\]", row.get("SNP", ""))
        if not m:
            continue
        snp_a, snp_b = m.group(1), m.group(2)

        # Genomic alleles on + strand
        strand = row[col_strand]
        gref = strand_normalize(row[col_ref], strand)
        galt = strand_normalize(row[col_alt], strand)

        # Complement of SNP alleles (for 'complement' decisions)
        snp_a_comp = complement(snp_a)
        snp_b_comp = complement(snp_b)

        chr_  = row[col_chr]
        pos   = row[col_pos]
        base_decision = decision.replace("indel_", "")

        # Match SNP alleles to genomic alleles to determine which is the SNP ref
        snp_ref = None
        geno_order = f"{gref},{galt}"

        if base_decision == "as_is":
            if snp_a == gref and snp_b == galt:
                snp_ref = snp_a
                snp_alleles = f"{snp_a},{snp_b}"
                geno_alleles = f"{gref},{galt}"
            elif snp_b == gref and snp_a == galt:
                snp_ref = snp_b
                snp_alleles = f"{snp_a},{snp_b}"
                geno_alleles = f"{galt},{gref}"
            else:
                errors += 1
                lines.append(f"Error\t{name}\t{decision}\t{snp_a},{snp_b}\t{gref},{galt}")
                continue
        elif base_decision == "complement":
            if snp_a_comp == gref and snp_b_comp == galt:
                snp_ref = snp_a
                snp_alleles = f"{snp_a},{snp_b}"
                geno_alleles = f"{gref},{galt}"
            elif snp_b_comp == gref and snp_a_comp == galt:
                snp_ref = snp_b
                snp_alleles = f"{snp_a},{snp_b}"
                geno_alleles = f"{galt},{gref}"
            else:
                errors += 1
                lines.append(f"Error\t{name}\t{decision}\t{snp_a},{snp_b}\t{gref},{galt}")
                continue
        else:
            errors += 1
            lines.append(f"Error\t{name}\t{decision}\t{snp_a},{snp_b}\t{gref},{galt}")
            continue

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

    # ── Filter 1: Unmapped (Strand == N/A) ──────────────────────────────────
    df_mapped = df[df[col_strand].isin(["+", "-"])].copy()
    qc_stats["After unmapped filter"] = len(df_mapped)
    print(f"[qc] After unmapped filter: {len(df_mapped):,} ({len(df) - len(df_mapped):,} removed)")

    # ── Filter 2: MAPQ filter ────────────────────────────────────────────────
    df_topseq_filtered = df_mapped[df_mapped["MAPQ_TopGenomicSeq"] >= args.mapq_topseq]
    df_mapq = apply_probe_mapq_filter(df_topseq_filtered, args.mapq_probe).copy()
    qc_stats["After MAPQ filter"] = len(df_mapq)
    print(f"[qc] After MAPQ filter (TopSeq>={args.mapq_topseq}): {len(df_mapq):,} ({len(df_mapped) - len(df_mapq):,} removed)")

    # ── Filter 2.5: CoordDelta filter ────────────────────────────────────────
    col_coord_delta  = f"CoordDelta_{assembly}"
    col_map_status   = f"MappingStatus_{assembly}"
    if args.coord_delta >= 0:
        if col_coord_delta not in df_mapq.columns or col_map_status not in df_mapq.columns:
            print(f"[qc] WARNING: --coord-delta requested but {col_coord_delta!r} or "
                  f"{col_map_status!r} column not found in input. Skipping filter. "
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
            is_topseq_only = df_mapq[col_map_status] == "topseq_only"
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

    # ── Allele usage decisions ───────────────────────────────────────────────
    print("[qc] Computing allele usage decisions...")
    decisions = {}
    decision_path = os.path.join(out_dir, "allele_usage_decision.txt")
    with open(decision_path, "w") as f:
        for _, row in df_coord.iterrows():
            d = allele_usage_decision(row, col_strand, assembly)
            decisions[row["Name"]] = d
            f.write(f"{row['Name']}\t{d}\n")
    print(f"[qc] Allele decisions written to: {decision_path}")

    # ── Generate VCF of all remapped positions ───────────────────────────────
    print("[qc] Generating VCF position template...")
    pos_vcf  = os.path.join(temp_dir, "_pos.vcf")
    ref_vcf  = os.path.join(temp_dir, "_ref.vcf")
    build_pos_vcf(df_coord, args.vcf_contigs, col_chr, col_pos, pos_vcf)

    print("[qc] Extracting reference alleles with bcftools...")
    ref_alleles = extract_ref_alleles(pos_vcf, args.reference, ref_vcf)

    # ── Strand-normalise remapped alleles (to + strand) ─────────────────────
    df_coord = df_coord.copy()
    df_coord["_gref"] = df_coord.apply(lambda r: strand_normalize(str(r[col_ref]), r[col_strand]), axis=1)
    df_coord["_galt"] = df_coord.apply(lambda r: strand_normalize(str(r[col_alt]), r[col_strand]), axis=1)
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
    matching = df_coord[
        (df_coord["_gref"] == df_coord["_genome_ref"]) &
        (df_coord["_galt"] != "-")  # remove remaining indels
    ].copy()
    qc_stats["After design conflict filter"] = len(matching)
    print(f"[qc] After design conflict filter: {len(matching):,} ({len(df_coord) - len(matching):,} removed)")

    # Write _matchingSNPs.vcf
    match_vcf = os.path.join(out_dir, "_matchingSNPs.vcf")
    with open(match_vcf, "w") as f:
        with open(args.vcf_contigs) as vc:
            f.write("##fileformat=VCFv4.3\n")
            f.write(vc.read())
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in matching.iterrows():
            f.write(f"{row[col_chr]}\t{row[col_pos]}\t{row['Name']}\t{row['_gref']}\t{row['_galt']}\t.\t.\t.\n")

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
    bim_path = os.path.join(out_dir, f"{prefix}_remapped_{assembly}.bim")
    bim = df_consistent[[col_chr, "Name", col_pos, "_gref", "_galt"]].copy()
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
    errors = build_final_map(df_consistent, decisions, assembly, map_path)
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
