"""
remap_manifest.py — Dual-alignment manifest remapper.

Takes an Illumina genotyping array manifest (CSV) and a target reference genome (FASTA),
and remaps each SNP marker to the new assembly using a context-aware dual-alignment strategy:

  1. Probe alignment (AlleleA_ProbeSeq, 50 bp) → high-precision coordinate
  2. TopGenomicSeq alignment (full context with [A/B]) → strand authority + Ref/Alt determination

Outputs the original manifest with these columns appended:
  Chr_{assembly}          Chromosome on the new assembly
  MapInfo_{assembly}      Base-pair position
  Strand_{assembly}       Alignment strand (+ / - / N/A)
  Ref_{assembly}          Reference allele
  Alt_{assembly}          Alternate allele
  MAPQ_TopGenomicSeq      Mapping quality of the winning TopGenomicSeq alignment
  MAPQ_Probe              Mapping quality of the selected probe alignment (NaN = no probe alignment, topseq_only)

Usage:
  python scripts/remap_manifest.py \\
      -i original_manifest.csv \\
      -r reference.fa \\
      -o output_remapped.csv \\
      -a equCab3 \\
      [--threads 4] \\
      [--temp-dir /tmp/remap]
"""

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass

import pandas as pd
import pysam

# IUPAC ambiguity codes (excluding N/n) → A; str.translate is O(n) C-level loop
_IUPAC_TO_A = str.maketrans("MRWSYKBDHVmrwsykbdhv", "A" * 20)

# Nucleotide complement lookup (single-base)
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Remap Illumina manifest probes to a new reference genome."
    )
    p.add_argument("-i", "--manifest", required=True, help="Input manifest CSV")
    p.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p.add_argument("-o", "--output", required=True, help="Output manifest CSV")
    p.add_argument(
        "-a", "--assembly",
        default="new_assembly",
        help="Assembly name used to label output columns and files (default: new_assembly)",
    )
    p.add_argument("--threads", type=int, default=4, help="Threads for minimap2 (default: 4)")
    p.add_argument(
        "--temp-dir",
        default=None,
        help="Directory for temporary FASTA/SAM files (default: same directory as output)",
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Skip minimap2 alignment if SAM files already exist in temp-dir",
    )
    return p.parse_args()

# ── MANIFEST PARSING ─────────────────────────────────────────────────────────

def locate_data_section(filename, start_marker="[Assay]", end_marker="[Controls]"):
    """Returns (skiprows, nrows) for the data section of an Illumina manifest."""
    header_line = 0
    footer_line = None
    with open(filename) as f:
        for i, line in enumerate(f):
            stripped = line.strip().rstrip(",")
            if stripped == start_marker:
                header_line = i + 1
            if stripped == end_marker:
                footer_line = i
                break
    nrows = (footer_line - header_line) if footer_line is not None else None
    return header_line, nrows

# ── SEQUENCE HELPERS ─────────────────────────────────────────────────────────

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence (supports upper and lower case)."""
    return seq.translate(_RC_TABLE)[::-1]


def extract_candidates(top_seq):
    """
    Parses TopGenomicSeq format 'PREFIX[A/B]SUFFIX' into (pre, alleleA, alleleB, post).
    Returns (None, None, None, None) if the format is not recognised.

    The deletion notation '-' (as in [-/CTCGTG] or [CTCGTG/-]) is normalised to ''
    (empty string) because '-' means 'no sequence', not a literal dash character.
    """
    m = re.search(r"(.*?)\[(.*?)/(.*?)\](.*)", top_seq or "")
    if m:
        a = "" if m.group(2) == "-" else m.group(2)
        b = "" if m.group(3) == "-" else m.group(3)
        return m.group(1), a, b, m.group(4)
    return None, None, None, None


def probe_topseq_orientation(probe_seq, topseq_a, topseq_b):
    """
    Determine the orientation of a probe relative to TopGenomicSeq by sequence comparison.

    Returns:
      "same"       — probe (as-is) is a substring of topseq_a or topseq_b
      "complement" — RC(probe) is a substring of topseq_a or topseq_b
      "unknown"    — neither matches
    """
    if probe_seq in topseq_a or probe_seq in topseq_b:
        return "same"
    rc = reverse_complement(probe_seq)
    if rc in topseq_a or rc in topseq_b:
        return "complement"
    return "unknown"


def compute_probe_strand_agreement(ilmn_strand, topseq_strand, probe_align_strand,
                                   probe_seq, topseq_a, topseq_b):
    """
    Compute the probe strand and whether it agrees with the IlmnStrand expectation.

    For IlmnStrand TOP/BOT: uses the probe alignment strand directly.
      TOP → probe strand expected to match TopSeq strand.
      BOT → probe strand expected to be opposite to TopSeq strand.

    For IlmnStrand PLUS/MINUS: uses sequence comparison (probe_topseq_orientation).
      PLUS → probe expected on '+' strand.
      MINUS → probe expected on '-' strand.

    Returns (probe_strand, agreement_as_expected) as strings:
      probe_strand          : '+', '-', or 'N/A'
      agreement_as_expected : 'True', 'False', or 'N/A'
    """
    ilmn = ilmn_strand.upper() if ilmn_strand else ""

    if ilmn in ("TOP", "BOT"):
        probe_strand = probe_align_strand
        if ilmn == "TOP":
            agreement = (probe_align_strand == topseq_strand)
        else:  # BOT
            agreement = (probe_align_strand != topseq_strand)
        return probe_strand, str(agreement)

    if ilmn in ("PLUS", "MINUS"):
        orientation = probe_topseq_orientation(probe_seq, topseq_a, topseq_b)
        if orientation == "unknown":
            return "N/A", "N/A"
        # Derive absolute probe strand from orientation + TopSeq strand
        if orientation == "same":
            probe_strand = topseq_strand
        else:  # complement
            probe_strand = "-" if topseq_strand == "+" else "+"
        expected_strand = "+" if ilmn == "PLUS" else "-"
        return probe_strand, str(probe_strand == expected_strand)

    return "N/A", "N/A"

# ── CIGAR UTILITIES ──────────────────────────────────────────────────────────

def _parse_cigar(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def cigar_ref_span(cigar):
    """Reference bases consumed by the alignment (for computing alignment end position)."""
    return sum(n for n, op in _parse_cigar(cigar) if op in "MDN=X")


def get_alignment_end(pos, cigar):
    """1-based end coordinate of an alignment (inclusive)."""
    return pos + cigar_ref_span(cigar) - 1


def parse_cigar_to_ref_pos(start_pos: int, cigar: str, query_index: int):
    """
    Maps a 0-based query index to a 1-based reference coordinate by walking the CIGAR.

    Returns (ref_pos, in_softclip):
      ref_pos     — 1-based reference position of the target query base
      in_softclip — True if the target falls inside a soft-clipped region
                    (ref_pos is the clip junction in that case, not exact)

    Used to derive the SNP coordinate from the TopGenomicSeq alignment independent
    of the probe alignment, for cross-validation only.  The target query index is:
      + strand: info["PreLen"]   (start of allele bracket in original query)
      − strand: info["PostLen"]  (start of allele bracket in RC query)

    Limitation: on the minus strand, this returns the reference position of the
    first base of RC(allele), which is the *rightmost* reference coordinate of the
    allele.  For single-nucleotide variants (allele_len == 1) this equals the SNP
    start and CoordDelta will be 0.  For multi-base indels on minus strand,
    CoordDelta may be inflated by up to allele_len − 1 bases relative to the
    probe-based coordinate; filter indel markers before using CoordDelta.
    """
    ops = _parse_cigar(cigar)
    curr_q = 0
    curr_r = start_pos
    for n, op in ops:
        if op in "M=X":
            if query_index < curr_q + n:
                return curr_r + (query_index - curr_q), False
            curr_q += n
            curr_r += n
        elif op == "I":
            if query_index < curr_q + n:
                return curr_r, False   # inside insertion: return junction
            curr_q += n
        elif op == "S":
            if query_index < curr_q + n:
                return curr_r, True    # inside soft clip: approximate
            curr_q += n
        elif op in "DN":
            curr_r += n
    return curr_r, False


def get_probe_coordinate(pos_start, cigar_str, strand, assay_type):
    """
    Calculates the variant start coordinate from a probe alignment.

    Infinium II: variant is the base AFTER the probe 3' end.
    Infinium I:  variant is the LAST base of the probe.

    On the minus strand, the probe's physical 3' end is at the alignment start (POS).
    On the plus strand, it is at alignment start + reference_span - 1.
    """
    ops = _parse_cigar(cigar_str)

    if strand == "+":
        ref_span = sum(n for n, op in ops if op in "MDN=XS")
        probe_end = pos_start + ref_span - 1
        return probe_end + 1 if assay_type == "II" else probe_end
    else:
        leading_s = ops[0][0] if ops and ops[0][1] == "S" else 0
        probe_end = pos_start - leading_s
        return probe_end - 1 if assay_type == "II" else probe_end


def calculate_overlap(s1, e1, s2, e2):
    """Overlap length between two closed intervals [s1,e1] and [s2,e2]."""
    return max(0, min(e1, e2) - max(s1, s2))


def compute_qcov(cigar: str) -> float:
    """
    Fraction of query bases covered by alignment matches (M/=/X operations).
    Insertions and soft clips count toward total query length but not toward coverage.
    Returns 0.0 for an empty or unrecognised CIGAR.
    """
    ops = _parse_cigar(cigar)
    aligned = sum(n for n, op in ops if op in "M=X")
    total   = sum(n for n, op in ops if op in "MIS=X")  # query-consuming ops (H excluded: not in SEQ)
    return aligned / total if total > 0 else 0.0


def compute_soft_clip_frac(cigar: str) -> float:
    """
    Fraction of query bases that are soft-clipped.
    Returns 0.0 if there are no soft clips.
    """
    ops = _parse_cigar(cigar)
    clipped = sum(n for n, op in ops if op == "S")
    total   = sum(n for n, op in ops if op in "MIS=X")
    return clipped / total if total > 0 else 0.0


# ── SAM PARSING ──────────────────────────────────────────────────────────────

def _get_nm(cols):
    for tag in cols[11:]:
        if tag.startswith("NM:i:"):
            return int(tag.split(":")[2])
    return 999


def _get_as(cols):
    for tag in cols[11:]:
        if tag.startswith("AS:i:"):
            return int(tag.split(":")[2])
    return -1


def parse_topseq_sam(sam_path):
    """
    Reads the TopGenomicSeq SAM and returns a dict:
      { snp_name: { 'A': [align_dict, ...], 'B': [align_dict, ...] } }

    Primary and secondary alignments are both retained to give the pair-selection
    algorithm the full set of candidate loci. Supplementary alignments (chimeric
    split-reads, FLAG & 2048) are still discarded.
    """
    results = {}
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@") or line.startswith("[M"):
                continue
            cols = line.split("\t")
            flag = int(cols[1])
            if flag & 4:    # unmapped
                continue
            if flag & 2048: # supplementary
                continue
            qname_full = cols[0]
            which = qname_full[-1]      # 'A' or 'B'
            qname = qname_full[:-2]     # strip '_A' or '_B'
            pos = int(cols[3])
            cigar = cols[5]
            entry = {
                "NM":     _get_nm(cols),
                "AS":     _get_as(cols),
                "Chr":    cols[2],
                "Pos":    pos,
                "Cigar":  cigar,
                "MAPQ":   int(cols[4]),
                "Strand": "-" if flag & 16 else "+",
                "End":    get_alignment_end(pos, cigar),
            }
            results.setdefault(qname, {"A": [], "B": []})[which].append(entry)
    return results


def parse_probe_sam(sam_path):
    """
    Reads the probe SAM and returns a dict:
      { snp_name: [ {Chr, Pos, Cigar, Strand, MAPQ, AS, NM, End}, ... ] }
    All mapped alignments are kept (primary + secondary, for overlap checking).
    """
    results = {}
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@") or line.startswith("[M"):
                continue
            cols = line.split("\t")
            flag = int(cols[1])
            if flag & 4:  # unmapped
                continue
            pos = int(cols[3])
            cigar = cols[5]
            entry = {
                "Chr": cols[2],
                "Pos": pos,
                "Cigar": cigar,
                "Strand": "-" if flag & 16 else "+",
                "MAPQ": int(cols[4]),
                "AS": _get_as(cols),
                "NM": _get_nm(cols),
                "End": get_alignment_end(pos, cigar),
            }
            results.setdefault(cols[0], []).append(entry)
    return results

# ── PAIR SELECTION ────────────────────────────────────────────────────────────

def is_placed_chromosome(name):
    """
    Returns True if *name* is a standard assembled chromosome rather than an
    unplaced scaffold.  Matches digits, X, Y, MT/M with an optional 'chr' prefix.
    No species-specific list is hardcoded.
    """
    return bool(re.match(r"^(chr)?(\d+|X|Y|MT|M)$", name, re.IGNORECASE))

def _make_competing_rows(pairs, reason):
    """
    Builds a list of row dicts for the ambiguous/scaffold CSV files.

    pairs  — list of (allele, topseq_align_dict, probe_align_dict)
    reason — 'position_tie' | 'NM_tie' | 'scaffold_resolved'
    """
    rows = []
    for rank, (allele, ts, pb) in enumerate(pairs, 1):
        rows.append({
            "AmbiguityReason": reason,
            "PairRank":        rank,
            "TopSeqAllele":    allele,
            "TopSeqChr":       ts["Chr"],
            "TopSeqPos":       ts["Pos"],
            "TopSeqStrand":    ts["Strand"],
            "TopSeqMAPQ":      ts["MAPQ"],
            "TopSeqNM":        ts["NM"],
            "ProbeChr":        pb["Chr"],
            "ProbePos":        pb["Pos"],
            "ProbeMAPQ":       pb["MAPQ"],
            "MinMAPQ":         min(ts["MAPQ"], pb["MAPQ"]),
        })
    return rows


def select_best_pair(topseq_aligns, probe_aligns):
    """
    Finds the best (TopGenomicSeq × probe) alignment pair for a single marker.

    Both TopGenomicSeq and probe may have multiple alignments (primary + secondary).
    A triple (topseq_allele, topseq_align, probe_align) is valid only when:
      1. chromosome matches
      2. overlap between the two alignment windows > 0

    No strand constraint is applied: probes designed for the bottom strand align
    to the opposite strand from TopGenomicSeq and are still valid.  The output
    marker strand is always taken from the TopGenomicSeq alignment; the probe
    strand is used only internally by get_probe_coordinate().

    Valid triples are ranked by min(topseq_MAPQ, probe_MAPQ) — weakest-link
    scoring — then by TopSeq AS score, then by combined MAPQ as a tertiary sort.

    Among triples tied at the top score, unique loci (TopSeq chr:pos) are
    extracted and resolved:
      • 1 unique locus            → ('winner', allele, topseq_align, probe_align)
      • placed chrom + scaffolds  → ('scaffold_resolved', allele, topseq_align,
                                      probe_align, competing_rows)
      • 2+ placed-chrom loci      → ('ambiguous', competing_rows)
      • no valid triples          → None  (unmapped)

    topseq_aligns : {'A': [align_dict, ...], 'B': [align_dict, ...]}
    probe_aligns  : [align_dict, ...]
    """
    valid = []
    for allele, ts_list in topseq_aligns.items():
        for ts in ts_list:
            for pb in probe_aligns:
                if ts["Chr"] != pb["Chr"]:
                    continue
                if calculate_overlap(ts["Pos"], ts["End"], pb["Pos"], pb["End"]) <= 0:
                    continue
                min_mapq = min(ts["MAPQ"], pb["MAPQ"])
                ts_as    = ts.get("AS", -1)
                sum_mapq = ts["MAPQ"] + pb["MAPQ"]
                valid.append((min_mapq, ts_as, sum_mapq, allele, ts, pb))

    if not valid:
        return None

    valid.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
    top_min_q, top_ts_as, top_sum_q = valid[0][0], valid[0][1], valid[0][2]
    top_pairs = [
        (allele, ts, pb)
        for min_q, t_as, sum_q, allele, ts, pb in valid
        if min_q == top_min_q and t_as == top_ts_as and sum_q == top_sum_q
    ]

    unique_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in top_pairs}

    if len(unique_loci) == 1:
        best = top_pairs[0]
        return ("winner", best[0], best[1], best[2])

    placed    = [(a, ts, pb) for a, ts, pb in top_pairs if     is_placed_chromosome(ts["Chr"])]
    scaffolds = [(a, ts, pb) for a, ts, pb in top_pairs if not is_placed_chromosome(ts["Chr"])]

    unique_placed_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in placed}

    if placed and scaffolds and len(unique_placed_loci) == 1:
        competing = _make_competing_rows(top_pairs, "scaffold_resolved")
        best = placed[0]
        return ("scaffold_resolved", best[0], best[1], best[2], competing)

    # NM tiebreaker: among placed-chromosome loci, pick the one with the lowest
    # per-locus minimum NM.  Only applies when at least one placed locus exists.
    if placed:
        locus_min_nm = {}
        for a, ts, pb in placed:
            locus = (ts["Chr"], ts["Pos"])
            if locus not in locus_min_nm or ts["NM"] < locus_min_nm[locus]:
                locus_min_nm[locus] = ts["NM"]
        best_nm = min(locus_min_nm.values())
        nm_winning_loci = {loc for loc, nm in locus_min_nm.items() if nm == best_nm}
        if len(nm_winning_loci) == 1:
            winning_chr, winning_pos = next(iter(nm_winning_loci))
            nm_pairs = [(a, ts, pb) for a, ts, pb in top_pairs
                        if ts["Chr"] == winning_chr and ts["Pos"] == winning_pos]
            competing = _make_competing_rows(top_pairs, "position_tie")
            best = nm_pairs[0]
            return ("nm_position_resolved", best[0], best[1], best[2], competing)

    competing = _make_competing_rows(top_pairs, "position_tie")
    return ("ambiguous", competing)


# ── TOPSEQ-ONLY RESCUE ───────────────────────────────────────────────────────

def _best_topseq_for_rescue(ts_aligns):
    """
    Pick the best mapped TopSeq alignment across both alleles for CIGAR-only rescue.
    Used when no valid (TopSeq × probe) triple exists but TopSeq itself aligned.

    Returns (allele, alignment_dict) for the highest-MAPQ then highest-AS alignment,
    or (None, None) if no mapped alignment exists.
    """
    candidates = []
    for allele, aligns in ts_aligns.items():
        for a in aligns:
            if a.get("Chr", "*") not in ("*", "0"):
                candidates.append((allele, a))
    if not candidates:
        return None, None
    candidates.sort(key=lambda x: (x[1]["MAPQ"], x[1].get("AS", -1)), reverse=True)
    return candidates[0]


# ── REF/ALT DETERMINATION ────────────────────────────────────────────────────

def determine_ref_alt(winning_allele, winning_ts, topseq_aligns, candidates_info):
    """
    Determines reference and alternate alleles by comparing NM (edit distance)
    between the two TopGenomicSeq allele sequences at the winning chromosome.

    The allele whose sequence is more similar to the reference genome (lower NM)
    is the reference allele.  If NM is equal the marker is declared ambiguous
    (returns None) rather than using an arbitrary tiebreak.

    winning_allele  : 'A' or 'B'
    winning_ts      : the TopGenomicSeq alignment dict that defined the locus
    topseq_aligns   : {'A': [align_dict, ...], 'B': [align_dict, ...]}
    candidates_info : {'AlleleA': <nucleotide>, 'AlleleB': <nucleotide>, ...}

    Returns (ref_char, alt_char) or None if NM is tied.
    """
    other_allele = "B" if winning_allele == "A" else "A"
    nm_winner = winning_ts["NM"]

    other_at_chr = [
        a for a in topseq_aligns.get(other_allele, [])
        if a["Chr"] == winning_ts["Chr"]
    ]
    nm_other = min((a["NM"] for a in other_at_chr), default=999)

    if nm_winner < nm_other:
        ref_allele = winning_allele
    elif nm_other < nm_winner:
        ref_allele = other_allele
    else:
        return None  # NM tie → ambiguous

    alt_allele = "B" if ref_allele == "A" else "A"
    ref_char = candidates_info["AlleleA"] if ref_allele == "A" else candidates_info["AlleleB"]
    alt_char = candidates_info["AlleleB"] if ref_allele == "A" else candidates_info["AlleleA"]
    return ref_char, alt_char


def resolve_ref_from_genome(fasta, chr_, var_pos, allele_a_char, allele_b_char, strand):
    """
    Fetch the reference base at var_pos and determine which allele is Ref.

    allele_a_char / allele_b_char are in alignment-strand orientation.
    For minus-strand markers these are complemented before comparison with the
    forward-strand genome base.

    fasta          : open pysam.FastaFile
    chr_           : chromosome name
    var_pos        : 1-based variant position
    allele_a_char  : single nucleotide for allele A (alignment strand)
    allele_b_char  : single nucleotide for allele B (alignment strand)
    strand         : '+' or '-'

    Returns (ref_char, alt_char) in alignment-strand orientation, or None.
    """
    try:
        ref_base = fasta.fetch(chr_, var_pos - 1, var_pos).upper()
    except (ValueError, KeyError):
        return None

    if strand == "-":
        cmp_a = _COMP.get(allele_a_char, allele_a_char)
        cmp_b = _COMP.get(allele_b_char, allele_b_char)
    else:
        cmp_a = allele_a_char
        cmp_b = allele_b_char

    if ref_base == cmp_a:
        return allele_a_char, allele_b_char
    if ref_base == cmp_b:
        return allele_b_char, allele_a_char
    return None


# ── DECISION COUNTERS ─────────────────────────────────────────────────────────

@dataclass
class DecisionCounters:
    """
    Accumulates per-marker decision counts throughout the pipeline and prints
    a structured summary table on completion.
    """
    total_loaded:            int = 0
    topseq_both:             int = 0
    topseq_one:              int = 0
    topseq_neither:          int = 0
    valid_pair_found:        int = 0
    no_valid_pair:           int = 0
    unique_position:         int = 0
    scaffold_resolved:       int = 0
    nm_position_resolved:    int = 0
    position_ambiguous:      int = 0
    ref_alt_clear:           int = 0
    ref_alt_ref_resolved:    int = 0
    ref_alt_nm_tie:          int = 0
    mapq_60:                 int = 0
    mapq_30_59:              int = 0
    mapq_1_29:               int = 0
    mapq_0:                  int = 0
    final_mapped:            int = 0
    final_ref_resolved:      int = 0
    final_scaffold_resolved:      int = 0
    final_nm_position_resolved:   int = 0
    final_unmapped:               int = 0
    final_ambiguous:         int = 0
    final_topseq_only:            int = 0
    topseq_rescue_failed_softclip: int = 0   # SNP target in soft-clipped region
    topseq_rescue_failed_refalt:   int = 0   # NM tie unresolvable by ref lookup
    ref_base_mismatch:             int = 0
    strand_agreement_unexpected:   int = 0   # StrandAgreementAsExpected == False

    def format_summary(self) -> str:
        W = 55
        lines = []

        def row(label, n):
            lines.append(f"  {label:<{W - 2}} {n:>8,}")

        lines.append("")
        lines.append("=== REMAPPING DECISION SUMMARY ===")
        lines.append(f"  {'Total markers loaded:':<{W - 2}} {self.total_loaded:>8,}")
        lines.append("")
        lines.append("TopGenomicSeq alignment:")
        row("Both A and B aligned:", self.topseq_both)
        row("Only one allele aligned:", self.topseq_one)
        row("Neither aligned \u2192 unmapped:", self.topseq_neither)
        lines.append("")
        lines.append("Pair filtering (chr + strand + overlap):")
        row("\u22651 valid pair:", self.valid_pair_found)
        row("No valid pair (total):", self.no_valid_pair)
        row("  \u251c rescued via TopSeq CIGAR (topseq_only):", self.final_topseq_only)
        row("  \u251c rescue failed \u2014 SNP in soft-clipped region:", self.topseq_rescue_failed_softclip)
        row("  \u2514 rescue failed \u2014 NM tie unresolvable:", self.topseq_rescue_failed_refalt)
        lines.append("")
        lines.append("Position resolution:")
        row("Unique best pair:", self.unique_position)
        row("Scaffold\u2192chromosome tiebreak resolved:", self.scaffold_resolved)
        row("NM position-resolved (nm_position_resolved_markers.csv):", self.nm_position_resolved)
        row("True tie \u2192 ambiguous:", self.position_ambiguous)
        lines.append("")
        lines.append("Ref/Alt (NM comparison):")
        row("Clear NM winner:", self.ref_alt_clear)
        row("NM tie resolved by ref lookup:", self.ref_alt_ref_resolved)
        row("NM tie \u2192 ambiguous (triallelic):", self.ref_alt_nm_tie)
        row("Genome ref base mismatches (see RefBaseMatch col):", self.ref_base_mismatch)
        row("Strand agreement unexpected (StrandAgreementAsExpected=False):", self.strand_agreement_unexpected)
        lines.append("")
        lines.append("MAPQ distribution (min of winning pair):")
        row("MAPQ = 60:", self.mapq_60)
        row("MAPQ 30\u201359:", self.mapq_30_59)
        row("MAPQ  1\u201329:", self.mapq_1_29)
        row("MAPQ = 0:", self.mapq_0)
        lines.append("")
        lines.append("Final (all in main output):")
        row("mapped:", self.final_mapped)
        row("ref-resolved (ref lookup):", self.final_ref_resolved)
        row("nm-position-resolved:", self.final_nm_position_resolved)
        row("scaffold-resolved (scaffold_resolved_markers.csv):", self.final_scaffold_resolved)
        row("topseq-only (CIGAR rescue, no probe):", self.final_topseq_only)
        row("unmapped  (Chr=0):", self.final_unmapped)
        row("ambiguous (Chr=0 + ambiguous_markers.csv):", self.final_ambiguous)
        lines.append("")
        return "\n".join(lines)

# ── CORE REMAPPING ───────────────────────────────────────────────────────────

def run_remapping(args):
    assembly = args.assembly
    col_chr      = f"Chr_{assembly}"
    col_pos      = f"MapInfo_{assembly}"
    col_strand   = f"Strand_{assembly}"
    col_ref      = f"Ref_{assembly}"
    col_alt      = f"Alt_{assembly}"
    col_mapq_ts  = "MAPQ_TopGenomicSeq"
    col_mapq_pb  = "MAPQ_Probe"

    # ── Output / temp paths ──────────────────────────────────────────────────
    out_dir = os.path.dirname(os.path.abspath(args.output)) or "."
    temp_dir = os.path.abspath(args.temp_dir) if args.temp_dir else out_dir
    os.makedirs(temp_dir, exist_ok=True)

    topseq_fasta = os.path.join(temp_dir, "temp_topseq.fasta")
    probe_fasta  = os.path.join(temp_dir, "temp_probes.fasta")
    topseq_sam   = os.path.join(temp_dir, "temp_topseq.sam")
    probe_sam    = os.path.join(temp_dir, "temp_probe.sam")

    # ── Load manifest ────────────────────────────────────────────────────────
    print(f"[remap] Reading manifest: {args.manifest}")
    skip_rows, nrows = locate_data_section(args.manifest)
    dtype_dict = {
        "AddressA_ID": "string", "AddressB_ID": "string",
        "GenomeBuild": "string", "Chr": "string", "MapInfo": "Int64",
        "SourceVersion": "string", "BeadSetID": "string",
    }
    df = pd.read_csv(
        args.manifest, dtype=dtype_dict,
        skiprows=skip_rows, nrows=nrows, low_memory=False,
    )
    df = df.dropna(subset=["Name", "AlleleA_ProbeSeq"])
    print(f"[remap] Loaded {len(df):,} markers.")

    # ── Determine assay type (Infinium I vs II) ──────────────────────────────
    assay_types = {
        row["Name"]: ("II" if pd.isna(row.get("AlleleB_ProbeSeq")) else "I")
        for _, row in df.iterrows()
    }

    # ── Generate FASTA files ─────────────────────────────────────────────────
    print("[remap] Generating FASTA files for alignment...")
    candidates_info = {}
    with open(topseq_fasta, "w") as ft, open(probe_fasta, "w") as fp:
        for _, row in df.iterrows():
            name = row["Name"]
            fp.write(f">{name}\n{row['AlleleA_ProbeSeq']}\n")

            raw_topseq = row.get("TopGenomicSeq", "")
            if raw_topseq:
                raw_topseq = raw_topseq.translate(_IUPAC_TO_A)
            pre, a, b, post = extract_candidates(raw_topseq)
            if pre is not None:
                seq_a = pre + (a if a != "-" else "") + post
                seq_b = pre + (b if b != "-" else "") + post
                ft.write(f">{name}_A\n{seq_a}\n")
                ft.write(f">{name}_B\n{seq_b}\n")
                candidates_info[name] = {
                    "PreLen": len(pre), "PostLen": len(post),
                    "AlleleA": a, "AlleleB": b,
                    "TopSeqA": seq_a, "TopSeqB": seq_b,
                }

    # ── Align with minimap2 ──────────────────────────────────────────────────
    t = args.threads
    sams_exist = os.path.exists(topseq_sam) and os.path.exists(probe_sam)
    if args.resume and sams_exist:
        print(f"[remap] --resume: reusing existing SAM files in {temp_dir}")
    else:
        if args.resume and not sams_exist:
            print(f"[remap] --resume: SAM files not found in {temp_dir}, running alignment")
        else:
            print("[remap] Aligning sequences with minimap2...")
        subprocess.check_call(
            # -N 5: retain up to 5 secondary alignments per query for multi-locus resolution
            f"minimap2 -ax sr -N 5 -t {t} {args.reference} {topseq_fasta} > {topseq_sam} 2>/dev/null",
            shell=True,
        )
        subprocess.check_call(
            # -N 5: allow up to 5 secondary alignments for overlap checking
            f"minimap2 -ax sr -N 5 -t {t} {args.reference} {probe_fasta} > {probe_sam} 2>/dev/null",
            shell=True,
        )

    # ── Step 4: Parse SAM files ──────────────────────────────────────────────
    print("[remap] Parsing alignments...")
    raw_topseq       = parse_topseq_sam(topseq_sam)
    probe_candidates = parse_probe_sam(probe_sam)

    # ── Step 5: Resolve coordinates for each marker ──────────────────────────
    print("[remap] Resolving coordinates...")
    col_status = f"MappingStatus_{assembly}"
    col_delta  = "DeltaScore_TopGenomicSeq"
    col_qcov   = "QueryCov_TopGenomicSeq"
    col_scfrac = "SoftClipFrac_TopGenomicSeq"
    col_cigar_coord  = f"Coord_TopSeqCIGAR_{assembly}"
    col_probe_coord  = f"CoordProbe_{assembly}"
    col_coord_delta  = f"CoordDelta_{assembly}"
    col_coord_source = f"CoordSource_{assembly}"
    col_ref_match    = f"RefBaseMatch_{assembly}"
    col_probe_strand = f"ProbeStrand_{assembly}"
    col_strand_agree = f"StrandAgreementAsExpected_{assembly}"
    new_cols = {
        col_chr: [], col_pos: [], col_strand: [],
        col_ref: [], col_alt: [],
        col_mapq_ts: [], col_mapq_pb: [],
        col_delta: [],
        col_qcov:    [],
        col_scfrac:  [],
        col_cigar_coord:  [],
        col_probe_coord:  [],
        col_coord_delta:  [],
        col_coord_source: [],
        col_ref_match: [],
        col_probe_strand: [],
        col_strand_agree: [],
        col_status: [],
    }
    ambiguous_rows = []
    scaffold_rows  = []
    nm_pos_rows    = []
    counters       = DecisionCounters()
    counters.total_loaded = len(df)

    fasta = pysam.FastaFile(args.reference)
    try:

        def _append_unmapped_cols(status):
            new_cols[col_chr].append("0")
            new_cols[col_pos].append(0)
            new_cols[col_strand].append("N/A")
            new_cols[col_ref].append("N")
            new_cols[col_alt].append("N")
            new_cols[col_mapq_ts].append(0)
            new_cols[col_mapq_pb].append(0)
            new_cols[col_delta].append(-1)
            new_cols[col_qcov].append(0.0)
            new_cols[col_scfrac].append(0.0)
            new_cols[col_cigar_coord].append(0)
            new_cols[col_probe_coord].append(0)
            new_cols[col_coord_delta].append(-1)
            new_cols[col_coord_source].append("N/A")
            new_cols[col_ref_match].append("N/A")
            new_cols[col_probe_strand].append("N/A")
            new_cols[col_strand_agree].append("N/A")
            new_cols[col_status].append(status)

        for _, row in df.iterrows():
            name      = row["Name"]
            ts_aligns = raw_topseq.get(name, {})
            pb_aligns = probe_candidates.get(name, [])
            info      = candidates_info.get(name)

            # ΔScore: AS_best − AS_2nd across all TopSeq alignments.
            # -1 = no alignments with valid AS, or only one alignment (uniquely placed).
            # 0  = two alignments with identical AS (ambiguous locus).
            # >0 = margin of superiority of best alignment over second-best.
            all_as = sorted(
                [a["AS"] for aligns in ts_aligns.values() for a in aligns if a.get("AS", -1) >= 0],
                reverse=True,
            )
            delta_score = (all_as[0] - all_as[1]) if len(all_as) >= 2 else -1

            has_a = bool(ts_aligns.get("A"))
            has_b = bool(ts_aligns.get("B"))
            if has_a and has_b:
                counters.topseq_both += 1
            elif has_a or has_b:
                counters.topseq_one += 1
            else:
                counters.topseq_neither += 1

            if not info or (not has_a and not has_b):
                _append_unmapped_cols("unmapped")
                counters.final_unmapped += 1
                continue

            result = select_best_pair(ts_aligns, pb_aligns)

            if result is None:
                counters.no_valid_pair += 1
                # Rescue: TopSeq aligned but no valid (TopSeq × probe) triple.
                # Derive coordinate purely from TopGenomicSeq CIGAR walk.
                best_allele, best_ts = _best_topseq_for_rescue(ts_aligns)
                if best_ts is None:
                    _append_unmapped_cols("unmapped")
                    counters.final_unmapped += 1
                    continue
                target_idx = info["PreLen"] if best_ts["Strand"] == "+" else info["PostLen"]
                cigar_coord, cigar_in_sc = parse_cigar_to_ref_pos(
                    best_ts["Pos"], best_ts["Cigar"], target_idx
                )
                if cigar_in_sc or cigar_coord == 0:
                    # TopSeq aligned but SNP target index falls in soft-clipped region;
                    # no reference coordinate can be derived from this alignment.
                    _append_unmapped_cols("unmapped")
                    counters.topseq_rescue_failed_softclip += 1
                    counters.final_unmapped += 1
                    continue
                ref_alt = determine_ref_alt(best_allele, best_ts, ts_aligns, info)
                if ref_alt is None:
                    ref_alt = resolve_ref_from_genome(
                        fasta, best_ts["Chr"], cigar_coord, info["AlleleA"], info["AlleleB"], best_ts["Strand"]
                    )
                if ref_alt is None:
                    _append_unmapped_cols("unmapped")
                    counters.topseq_rescue_failed_refalt += 1
                    counters.final_unmapped += 1
                    continue
                ref_char, alt_char = ref_alt
                if len(ref_char) == 1 and len(alt_char) == 1:
                    try:
                        genome_base = fasta.fetch(
                            best_ts["Chr"], cigar_coord - 1, cigar_coord
                        ).upper()
                    except (ValueError, KeyError):
                        genome_base = None
                    ref_char_fwd = (
                        _COMP.get(ref_char, ref_char)
                        if best_ts["Strand"] == "-" else ref_char
                    )
                    if genome_base is None or genome_base == "":
                        ref_base_match_str = "False"
                        counters.ref_base_mismatch += 1
                    elif genome_base == ref_char_fwd:
                        ref_base_match_str = "True"
                    else:
                        ref_base_match_str = "False"
                        counters.ref_base_mismatch += 1
                else:
                    ref_base_match_str = "N/A"
                new_cols[col_chr].append(best_ts["Chr"])
                new_cols[col_pos].append(cigar_coord)
                new_cols[col_strand].append(best_ts["Strand"])
                new_cols[col_ref].append(ref_char)
                new_cols[col_alt].append(alt_char)
                new_cols[col_mapq_ts].append(best_ts["MAPQ"])
                new_cols[col_mapq_pb].append(float('nan'))  # no probe alignment
                new_cols[col_delta].append(delta_score)
                new_cols[col_qcov].append(compute_qcov(best_ts["Cigar"]))
                new_cols[col_scfrac].append(compute_soft_clip_frac(best_ts["Cigar"]))
                new_cols[col_cigar_coord].append(cigar_coord)
                new_cols[col_probe_coord].append(0)    # no probe alignment for this marker
                new_cols[col_coord_delta].append(-1)   # no probe coord to compare
                new_cols[col_coord_source].append("cigar")
                new_cols[col_ref_match].append(ref_base_match_str)
                new_cols[col_probe_strand].append("N/A")   # no probe alignment
                new_cols[col_strand_agree].append("N/A")   # no probe alignment
                new_cols[col_status].append("topseq_only")
                counters.final_topseq_only += 1
                continue

            counters.valid_pair_found += 1

            if result[0] == "ambiguous":
                _, competing = result
                ambiguous_rows.extend({"Name": name, **r} for r in competing)
                _append_unmapped_cols("ambiguous")
                counters.position_ambiguous += 1
                counters.final_ambiguous    += 1
                continue

            if result[0] == "scaffold_resolved":
                _, winning_allele, winning_ts, winning_pb, competing = result
                scaffold_rows.extend({"Name": name, **r} for r in competing)
                counters.scaffold_resolved += 1
            elif result[0] == "nm_position_resolved":
                _, winning_allele, winning_ts, winning_pb, competing = result
                nm_pos_rows.extend({"Name": name, **r} for r in competing)
                counters.nm_position_resolved += 1
            else:  # 'winner'
                _, winning_allele, winning_ts, winning_pb = result
                counters.unique_position += 1

            ref_alt = determine_ref_alt(winning_allele, winning_ts, ts_aligns, info)

            assay = assay_types.get(name, "II")
            c_pos = get_probe_coordinate(
                winning_pb["Pos"], winning_pb["Cigar"], winning_pb["Strand"], assay
            )

            if ref_alt is None:
                # Attempt to break NM tie by reference genome lookup at the variant position
                ref_alt = resolve_ref_from_genome(
                    fasta, winning_ts["Chr"], c_pos, info["AlleleA"], info["AlleleB"], winning_ts["Strand"]
                )
                if ref_alt is None:
                    # True triallelic — record both allele alignments for transparency
                    # (winning_pb reused for competing pair: no probe alignment for other allele)
                    other_allele = "B" if winning_allele == "A" else "A"
                    other_at_chr = [a for a in ts_aligns.get(other_allele, [])
                                    if a["Chr"] == winning_ts["Chr"]]
                    pairs = [(winning_allele, winning_ts, winning_pb)]
                    if other_at_chr:
                        other_ts = min(other_at_chr, key=lambda a: a["NM"])
                        pairs.append((other_allele, other_ts, winning_pb))
                    ambiguous_rows.extend({"Name": name, **r}
                                          for r in _make_competing_rows(pairs, "NM_tie"))
                    _append_unmapped_cols("ambiguous")
                    counters.ref_alt_nm_tie  += 1
                    counters.final_ambiguous += 1
                    continue
                counters.ref_alt_ref_resolved += 1
                status = "ref_resolved"
            else:
                counters.ref_alt_clear += 1
                if len(ref_alt[0]) > len(ref_alt[1]) and winning_ts["Strand"] == "-":
                    c_pos -= len(ref_alt[0]) - len(ref_alt[1])
                if result[0] in ("scaffold_resolved", "nm_position_resolved"):
                    status = result[0]
                else:
                    status = "mapped"

            ref_char, alt_char = ref_alt

            # Validate ref allele against genome reference base at the variant position.
            # For multi-character alleles (indels), skip validation — a single-base lookup
            # is not meaningful.
            #
            # ref_char is in the alignment strand (TopGenomicSeq orientation).  To compare
            # it against the forward-strand genome base we strand-normalise it first:
            # complement on minus-strand alignments, identity on plus.  This mirrors the
            # strand_normalize logic in qc_filter.py and produces the correct ~184 mismatches
            # that correspond to the markers qc_filter.py removes at the design-conflict step.
            if len(ref_char) == 1 and len(alt_char) == 1:
                try:
                    genome_base = fasta.fetch(
                        winning_ts["Chr"], c_pos - 1, c_pos
                    ).upper()
                except (ValueError, KeyError):
                    genome_base = None
                ref_char_fwd = (
                    _COMP.get(ref_char, ref_char)
                    if winning_ts["Strand"] == "-"
                    else ref_char
                )
                if genome_base is None or genome_base == "":
                    ref_base_match_str = "False"
                    counters.ref_base_mismatch += 1
                elif genome_base == ref_char_fwd:
                    ref_base_match_str = "True"
                else:
                    # Forward-strand genome base does not match strand-normalised ref.
                    ref_base_match_str = "False"
                    counters.ref_base_mismatch += 1
            else:
                ref_base_match_str = "N/A"  # indel: single-base lookup not applicable

            min_mapq = min(winning_ts["MAPQ"], winning_pb["MAPQ"])
            if min_mapq == 60:
                counters.mapq_60    += 1
            elif min_mapq >= 30:
                counters.mapq_30_59 += 1
            elif min_mapq >= 1:
                counters.mapq_1_29  += 1
            else:
                counters.mapq_0     += 1

            # CIGAR-based coordinate: cross-validator and coordinate override.
            # CoordDelta stores |probe_coord - cigar_coord| for diagnostics.
            # When delta >= 2, CIGAR coordinate is used as final position
            # (benchmark shows it outperforms probe coord at delta >= 2).
            # For indel markers (one allele is empty string), CIGAR coordinate
            # is always used as the final position when available: probe-based
            # coordinates are unreliable for multi-base alleles.
            target_idx = info["PreLen"] if winning_ts["Strand"] == "+" else info["PostLen"]
            cigar_coord, cigar_in_sc = parse_cigar_to_ref_pos(
                winning_ts["Pos"], winning_ts["Cigar"], target_idx
            )
            is_indel = len(ref_char) == 0 or len(alt_char) == 0
            if cigar_in_sc:
                cigar_out       = 0
                coord_delta_val = -1
                final_pos       = c_pos
                coord_source    = "probe"
            else:
                coord_delta_val = abs(c_pos - cigar_coord)
                cigar_out       = cigar_coord
                if is_indel or coord_delta_val >= 2:
                    final_pos    = cigar_coord
                    coord_source = "cigar"
                else:
                    final_pos    = c_pos
                    coord_source = "probe"

            ilmn_strand = str(row.get("IlmnStrand", "")).strip().upper()
            pb_strand, sa_expected = compute_probe_strand_agreement(
                ilmn_strand=ilmn_strand,
                topseq_strand=winning_ts["Strand"],
                probe_align_strand=winning_pb["Strand"],
                probe_seq=str(row.get("AlleleA_ProbeSeq", "")),
                topseq_a=info.get("TopSeqA", ""),
                topseq_b=info.get("TopSeqB", ""),
            )
            if sa_expected == "False":
                counters.strand_agreement_unexpected += 1

            new_cols[col_chr].append(winning_ts["Chr"])
            new_cols[col_pos].append(final_pos)
            new_cols[col_strand].append(winning_ts["Strand"])
            new_cols[col_ref].append(ref_char)
            new_cols[col_alt].append(alt_char)
            new_cols[col_mapq_ts].append(winning_ts["MAPQ"])
            new_cols[col_mapq_pb].append(winning_pb["MAPQ"])
            new_cols[col_delta].append(delta_score)
            new_cols[col_qcov].append(compute_qcov(winning_ts["Cigar"]))
            new_cols[col_scfrac].append(compute_soft_clip_frac(winning_ts["Cigar"]))
            new_cols[col_cigar_coord].append(cigar_out)
            new_cols[col_probe_coord].append(c_pos)   # raw probe-derived coord, before override
            new_cols[col_coord_delta].append(coord_delta_val)
            new_cols[col_coord_source].append(coord_source)
            new_cols[col_ref_match].append(ref_base_match_str)
            new_cols[col_probe_strand].append(pb_strand)
            new_cols[col_strand_agree].append(sa_expected)
            new_cols[col_status].append(status)

            if status == "scaffold_resolved":
                counters.final_scaffold_resolved += 1
            elif status == "nm_position_resolved":
                counters.final_nm_position_resolved += 1
            elif status == "ref_resolved":
                counters.final_ref_resolved += 1
            else:
                counters.final_mapped += 1

        # ── Step 6: Write outputs ────────────────────────────────────────────────
        for col, vals in new_cols.items():
            df[col] = vals

        print(f"[remap] Writing output: {args.output}")
        df.to_csv(args.output, index=False)

        if ambiguous_rows:
            amb_path = os.path.join(out_dir, "ambiguous_markers.csv")
            pd.DataFrame(ambiguous_rows).to_csv(amb_path, index=False)
            print(f"[remap] Ambiguous markers written to: {amb_path}")

        if scaffold_rows:
            scaf_path = os.path.join(out_dir, "scaffold_resolved_markers.csv")
            pd.DataFrame(scaffold_rows).to_csv(scaf_path, index=False)
            print(f"[remap] Scaffold-resolved markers written to: {scaf_path}")

        if nm_pos_rows:
            nm_pos_path = os.path.join(out_dir, "nm_position_resolved_markers.csv")
            pd.DataFrame(nm_pos_rows).to_csv(nm_pos_path, index=False)
            print(f"[remap] NM-position-resolved markers written to: {nm_pos_path}")

        print(f"[remap] Temp files in: {temp_dir}")

        # ── Step 7: Write decision summary to stdout and remapping_Report.txt ──
        summary = counters.format_summary()
        print(summary)
        report_path = os.path.join(out_dir, "remapping_Report.txt")
        with open(report_path, "w") as f:
            f.write(summary + "\n")
        print(f"[remap] Remapping report: {report_path}")

    finally:
        fasta.close()


if __name__ == "__main__":
    run_remapping(parse_args())
