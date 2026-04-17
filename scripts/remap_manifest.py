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
from collections import defaultdict
from dataclasses import dataclass

import pandas as pd
import pysam

# IUPAC ambiguity codes (excluding N/n) → A; str.translate is O(n) C-level loop
_IUPAC_TO_A = str.maketrans("MRWSYKBDHVmrwsykbdhv", "A" * 20)

from _strand_utils import strand_normalize

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

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence (upper or lower case).

    Thin wrapper over ``_strand_utils.strand_normalize(seq, "-")`` so that the
    reverse-complement table lives in exactly one place in the codebase.
    """
    return strand_normalize(seq, "-")


def _kmer_orientation(probe_seq, topseq_a, k=21):
    """
    Sequence-based fallback when substring matching cannot decide orientation.

    Builds k-mer sets for probe_seq, reverse_complement(probe_seq), and topseq_a,
    and returns whichever probe orientation shares more k-mers with topseq_a.
    Uses topseq_a only — at k=21 topseq_b differs only at the variant base and
    would not change the decision.

    Always returns "same" or "complement"; ties and degenerate inputs resolve
    to "same".
    """
    if not probe_seq or not topseq_a or len(probe_seq) < k or len(topseq_a) < k:
        return "same"
    rc = reverse_complement(probe_seq)
    topseq_kmers = {topseq_a[i:i + k] for i in range(len(topseq_a) - k + 1)}
    fwd_kmers    = {probe_seq[i:i + k] for i in range(len(probe_seq) - k + 1)}
    rc_kmers     = {rc[i:i + k]       for i in range(len(rc) - k + 1)}
    fwd_hits = len(fwd_kmers & topseq_kmers)
    rc_hits  = len(rc_kmers  & topseq_kmers)
    return "same" if fwd_hits >= rc_hits else "complement"


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

    Fast path — substring presence in topseq_a or topseq_b:
      "same"       — probe (as-is) is a substring
      "complement" — RC(probe) is a substring
    Fallback — k-mer overlap against topseq_a (k=21) via _kmer_orientation.

    Always returns "same" or "complement" — never "unknown".
    """
    if probe_seq in topseq_a or probe_seq in topseq_b:
        return "same"
    rc = reverse_complement(probe_seq)
    if rc in topseq_a or rc in topseq_b:
        return "complement"
    return _kmer_orientation(probe_seq, topseq_a)


def compute_probe_strand_agreement(topseq_strand, probe_align_strand,
                                   probe_seq, topseq_a, topseq_b):
    """
    Compute probe strand and whether it agrees with the sequence-derived expectation.

    The expected probe-vs-TopSeq strand relationship is derived purely from
    sequence comparison via probe_topseq_orientation (substring fast path +
    21-mer fallback against topseq_a). IlmnStrand is not consulted.

      orientation == "same"       → expected_probe_strand = topseq_strand
      orientation == "complement" → expected_probe_strand = flip(topseq_strand)

    Returns (probe_strand, agreement_as_expected):
      probe_strand          : '+' or '-' — pass-through of probe_align_strand
      agreement_as_expected : 'True' or 'False' — never 'N/A'
    """
    orientation = probe_topseq_orientation(probe_seq, topseq_a, topseq_b)
    if orientation == "same":
        expected_strand = topseq_strand
    else:  # "complement"
        expected_strand = "-" if topseq_strand == "+" else "+"
    return probe_align_strand, str(probe_align_strand == expected_strand)

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


# ── ALIGNMENT STATUS ─────────────────────────────────────────────────────────

def compute_alignment_status(ts_aligns, probe_aligns):
    """
    Raw alignment census — which sources produced at least one mapped hit.
    Called before any filtering or decision logic; 'aligned' means any mapped
    hit exists (any MAPQ, any chromosome).

    Returns one of: 'gp1', 'gp2', 'gp3', 'gp4', 'gp5', 'unmapped'.

    gp1: both TopSeq alleles + probe aligned
    gp2: exactly one TopSeq allele + probe aligned
    gp3: both TopSeq alleles, no probe
    gp4: exactly one TopSeq allele, no probe
    gp5: probe only (no TopSeq)
    unmapped: nothing aligned
    """
    has_a     = bool(ts_aligns.get("A"))
    has_b     = bool(ts_aligns.get("B"))
    has_probe = bool(probe_aligns)
    both_ts   = has_a and has_b
    one_ts    = has_a ^ has_b  # XOR: exactly one

    if both_ts and has_probe:
        return "gp1"
    if one_ts and has_probe:
        return "gp2"
    if both_ts and not has_probe:
        return "gp3"
    if one_ts and not has_probe:
        return "gp4"
    if has_probe and not has_a and not has_b:
        return "gp5"
    return "unmapped"


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
    Build row dicts for the *_unresolved_markers.csv / scaffold_resolved CSVs.

    pairs  — list of (allele, topseq_align_dict, probe_align_dict)
    reason — 'position_tie' | 'NM_tie' | 'AS_tie' | 'dAS_tie' |
             'CoordDelta_tie' | 'scaffold_resolved'
    """
    rows = []
    for rank, (allele, ts, pb) in enumerate(pairs, 1):
        rows.append({
            "UnresolvedReason": reason,
            "PairRank":         rank,
            "TopSeqAllele":     allele,
            "TopSeqChr":        ts["Chr"],
            "TopSeqPos":        ts["Pos"],
            "TopSeqStrand":     ts["Strand"],
            "TopSeqMAPQ":       ts["MAPQ"],
            "TopSeqNM":         ts["NM"],
            "ProbeChr":         pb["Chr"],
            "ProbePos":         pb["Pos"],
            "ProbeMAPQ":        pb["MAPQ"],
            "MinMAPQ":          min(ts["MAPQ"], pb["MAPQ"]),
        })
    return rows


def build_valid_triples(ts_aligns, probe_aligns,
                        probe_seq, topseq_a, topseq_b):
    """
    Enumerate valid (TopSeq_allele, ts_align, probe_align) triples.

    Validity requires:
      1. ts and probe on the same chromosome.
      2. Strand agreement (sequence-derived): compute_probe_strand_agreement
         must return 'True'. Expected probe strand is derived from
         probe_topseq_orientation + topseq_strand (21-mer fallback guarantees
         a decision — no 'N/A' is returned).
      3. Among strand-valid probes on the same chromosome, keep only the one
         with the highest overlap with ts (overlap-max selection).
      4. overlap(ts, best_probe) > 0.

    Returns a list of (allele, ts_align, probe_align) tuples.
    """
    triples = []
    for allele, ts_list in ts_aligns.items():
        for ts in ts_list:
            # Collect probes on the same chromosome
            same_chr = [pb for pb in probe_aligns if pb["Chr"] == ts["Chr"]]
            if not same_chr:
                continue
            # Strand-filter: keep probes where sequence-derived agreement is True
            strand_valid = []
            for pb in same_chr:
                _, agreement = compute_probe_strand_agreement(
                    topseq_strand=ts["Strand"],
                    probe_align_strand=pb["Strand"],
                    probe_seq=probe_seq,
                    topseq_a=topseq_a,
                    topseq_b=topseq_b,
                )
                if agreement == "True":
                    strand_valid.append(pb)
            if not strand_valid:
                continue
            # Overlap-max: keep the single strand-valid probe with highest overlap
            best_pb = max(
                strand_valid,
                key=lambda pb: calculate_overlap(ts["Pos"], ts["End"],
                                                  pb["Pos"], pb["End"])
            )
            if calculate_overlap(ts["Pos"], ts["End"],
                                  best_pb["Pos"], best_pb["End"]) <= 0:
                continue
            triples.append((allele, ts, best_pb))
    return triples


def _rank_single_aligns(candidates):
    """
    Rank a list of (label, align_dict) by AS → ΔAS → NM → scaffold_resolved.

    candidates : list of (label, align_dict) — only mapped alignments.
    Returns (label, align_dict, tie_status, competing_pool):
      - competing_pool is the list of (label, align_dict) pairs still tied
        after the final ranking step when tie_status == 'locus_unresolved';
        otherwise None.

    tie_status values: 'unique', 'AS_resolved', 'dAS_resolved', 'NM_resolved',
                       'scaffold_resolved', 'locus_unresolved', 'N/A' (no mapped aligns).
    """
    mapped = [(lbl, a) for lbl, a in candidates
              if a.get("Chr", "*") not in ("*", "0")]
    if not mapped:
        return None, None, "N/A", None

    # Unique locus check
    unique_loci = {(a["Chr"], a["Pos"]) for _, a in mapped}
    if len(unique_loci) == 1:
        return mapped[0][0], mapped[0][1], "unique", None

    all_as = [a.get("AS", -1) for _, a in mapped]

    # Step 1: AS — keep highest
    top_as = max(all_as)
    top = [(lbl, a) for lbl, a in mapped if a.get("AS", -1) == top_as]
    top_loci = {(a["Chr"], a["Pos"]) for _, a in top}
    if len(top_loci) == 1:
        return top[0][0], top[0][1], "AS_resolved", None

    # Step 2: ΔAS — for each surviving alignment, compute AS_this minus
    # the best AS of any OTHER alignment in the full mapped pool.
    def _das(a):
        other = max(
            (x.get("AS", -1) for _, x in mapped
             if (x["Chr"], x["Pos"]) != (a["Chr"], a["Pos"])),
            default=None,
        )
        return a.get("AS", -1) - other if other is not None else a.get("AS", -1)

    top_das = max(_das(a) for _, a in top)
    top = [(lbl, a) for lbl, a in top if _das(a) == top_das]
    top_loci = {(a["Chr"], a["Pos"]) for _, a in top}
    if len(top_loci) == 1:
        return top[0][0], top[0][1], "dAS_resolved", None

    # Step 3: NM — keep lowest
    min_nm = min(a.get("NM", 999) for _, a in top)
    top = [(lbl, a) for lbl, a in top if a.get("NM", 999) == min_nm]
    top_loci = {(a["Chr"], a["Pos"]) for _, a in top}
    if len(top_loci) == 1:
        return top[0][0], top[0][1], "NM_resolved", None

    # Step 4: scaffold_resolved — prefer placed chromosome over scaffold
    placed = [(lbl, a) for lbl, a in top if is_placed_chromosome(a["Chr"])]
    placed_loci = {(a["Chr"], a["Pos"]) for _, a in placed}
    if placed and len(placed_loci) == 1:
        return placed[0][0], placed[0][1], "scaffold_resolved", None

    return None, None, "locus_unresolved", top


def best_topseq_rescue(ts_aligns):
    """
    Pick the best TopSeq alignment across both alleles when no valid triple exists.

    Uses _rank_single_aligns (AS → ΔAS → NM → scaffold_resolved → locus_unresolved).
    Returns (allele, align_dict, tie_status, competing_pool):
      competing_pool is populated (list of (allele, align_dict)) only when
      tie_status == "locus_unresolved"; otherwise None.
    """
    candidates = [
        (allele, a)
        for allele, aligns in ts_aligns.items()
        for a in aligns
        if a.get("Chr", "*") not in ("*", "0")
    ]
    return _rank_single_aligns(candidates)


def best_probe_rescue(probe_aligns):
    """
    Pick the best probe alignment when TopSeq did not align at all.

    Uses _rank_single_aligns (AS → ΔAS → NM → scaffold_resolved → locus_unresolved).
    No strand filtering — without a TopSeq strand anchor, expected strand
    cannot be determined for TOP/BOT markers.
    Returns (align_dict, tie_status, competing_pool):
      competing_pool is populated only when tie_status == "locus_unresolved".
    """
    candidates = [
        ("probe", pb)
        for pb in probe_aligns
        if pb.get("Chr", "*") not in ("*", "0")
    ]
    _, align, tie, pool = _rank_single_aligns(candidates)
    return align, tie, pool


def rank_and_resolve(triples, all_ts_aligns, all_pb_aligns, info, assay_type):
    """
    Rank valid (allele, ts_align, probe_align) triples and resolve to a winner.

    Ranking steps (applied in order; each step only runs if previous left a tie):
      0. Unique locus: all triples point to same chr:pos → 'unique'
      1. AS sum: ts.AS + pb.AS, higher wins → 'AS_resolved'
      2. ΔAS sum: ΔAS_ts + ΔAS_pb, higher wins → 'dAS_resolved'
         ΔAS for each align = AS_this - max(AS of other aligns not at this locus)
      3. NM sum: ts.NM + pb.NM, lower wins → 'NM_resolved'
      4. CoordDelta: |probe_coord - CIGAR_coord|, lower wins → 'CoordDelta_resolved'
      5. Scaffold resolved: placed chr over scaffold → 'scaffold_resolved'
      6. Unresolved → 'locus_unresolved'

    Returns one of:
      ('unique'|'AS_resolved'|'dAS_resolved'|'NM_resolved'|
       'CoordDelta_resolved'|'scaffold_resolved', allele, ts, pb [,competing])
      ('locus_unresolved', competing)
    """
    all_ts_flat = [a for aligns in all_ts_aligns.values() for a in aligns]
    all_pb_flat = list(all_pb_aligns)

    # Step 0: unique locus check
    unique_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in triples}
    if len(unique_loci) == 1:
        allele, ts, pb = triples[0]
        return ("unique", allele, ts, pb)

    # Step 1: AS sum
    def _as_sum(triple):
        _, ts, pb = triple
        return ts.get("AS", -1) + pb.get("AS", -1)

    top_as = max(_as_sum(t) for t in triples)
    top = [t for t in triples if _as_sum(t) == top_as]
    top_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in top}
    if len(top_loci) == 1:
        competing = _make_competing_rows(top, "AS_tie")
        allele, ts, pb = top[0]
        return ("AS_resolved", allele, ts, pb, competing)

    # Step 2: ΔAS sum
    def _das_sum(triple):
        _, ts, pb = triple
        other_ts = max(
            (a.get("AS", -1) for a in all_ts_flat
             if (a["Chr"], a["Pos"]) != (ts["Chr"], ts["Pos"])),
            default=None,
        )
        das_ts = ts.get("AS", -1) - other_ts if other_ts is not None else ts.get("AS", -1)
        other_pb = max(
            (a.get("AS", -1) for a in all_pb_flat
             if (a["Chr"], a["Pos"]) != (pb["Chr"], pb["Pos"])),
            default=None,
        )
        das_pb = pb.get("AS", -1) - other_pb if other_pb is not None else pb.get("AS", -1)
        return das_ts + das_pb

    top_das = max(_das_sum(t) for t in top)
    top = [t for t in top if _das_sum(t) == top_das]
    top_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in top}
    if len(top_loci) == 1:
        competing = _make_competing_rows(top, "dAS_tie")
        allele, ts, pb = top[0]
        return ("dAS_resolved", allele, ts, pb, competing)

    # Step 3: NM sum
    def _nm_sum(triple):
        _, ts, pb = triple
        return ts.get("NM", 999) + pb.get("NM", 999)

    min_nm = min(_nm_sum(t) for t in top)
    top = [t for t in top if _nm_sum(t) == min_nm]
    top_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in top}
    if len(top_loci) == 1:
        competing = _make_competing_rows(top, "NM_tie")
        allele, ts, pb = top[0]
        return ("NM_resolved", allele, ts, pb, competing)

    # Step 4: CoordDelta — compute for all surviving triples
    def _coord_delta(triple):
        _, ts, pb = triple
        target_idx = info["PreLen"] if ts["Strand"] == "+" else info["PostLen"]
        c_pos = get_probe_coordinate(pb["Pos"], pb["Cigar"], pb["Strand"], assay_type)
        cigar_coord, in_sc = parse_cigar_to_ref_pos(ts["Pos"], ts["Cigar"], target_idx)
        if in_sc or cigar_coord == 0:
            return float("inf")  # unavailable → treat as worst
        return abs(c_pos - cigar_coord)

    min_delta = min(_coord_delta(t) for t in top)
    if min_delta < float("inf"):
        top_cd = [t for t in top if _coord_delta(t) == min_delta]
        top_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in top_cd}
        if len(top_loci) == 1:
            competing = _make_competing_rows(top, "CoordDelta_tie")
            allele, ts, pb = top_cd[0]
            return ("CoordDelta_resolved", allele, ts, pb, competing)
        top = top_cd

    # Step 5: scaffold_resolved
    placed    = [(a, ts, pb) for a, ts, pb in top if     is_placed_chromosome(ts["Chr"])]
    scaffolds = [(a, ts, pb) for a, ts, pb in top if not is_placed_chromosome(ts["Chr"])]
    placed_loci = {(ts["Chr"], ts["Pos"]) for _, ts, _ in placed}
    if placed and scaffolds and len(placed_loci) == 1:
        competing = _make_competing_rows(top, "scaffold_resolved")
        allele, ts, pb = placed[0]
        return ("scaffold_resolved", allele, ts, pb, competing)

    # Step 6: locus unresolved
    competing = _make_competing_rows(top, "position_tie")
    return ("locus_unresolved", competing)


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

    cmp_a = strand_normalize(allele_a_char, strand)
    cmp_b = strand_normalize(allele_b_char, strand)

    if ref_base == cmp_a:
        return allele_a_char, allele_b_char
    if ref_base == cmp_b:
        return allele_b_char, allele_a_char
    return None


def _refine_deletion_pos(fasta, chrom, initial_pos, gref_fwd, max_offset=10):
    """
    Search within ±max_offset bases of initial_pos for the exact start of gref_fwd
    in the genome. Scans outward one base at a time: 0, +1, -1, +2, -2, ...
    Returns the matching 1-based position, or None if not found within the window.
    """
    for offset in range(0, max_offset + 1):
        for delta in ([0] if offset == 0 else [offset, -offset]):
            pos = initial_pos + delta
            if pos < 1:
                continue
            try:
                if fasta.fetch(chrom, pos - 1, pos - 1 + len(gref_fwd)).upper() \
                        == gref_fwd.upper():
                    return pos
            except (ValueError, KeyError):
                pass
    return None


def determine_ref_alt_v2(winning_allele, winning_ts, ts_aligns,
                          candidates_info, fasta, chr_, final_pos, strand):
    """
    Determine Ref/Alt alleles for a mapped marker.

    Called after coordinate computation so final_pos = MapInfo (post-CoordDelta).

    For SNPs (both alleles are single nucleotides):
      - Genome lookup (primary): resolve_ref_from_genome → strand-aware.
      - NM comparison (parallel): determine_ref_alt when TopSeq aligned.
      - agreement_str: 'NM_match'|'NM_unmatch'|'NM_tied'|'NM_N/A'|'NM_only'|'refalt_unresolved'

    For indels (at least one allele is multi-base or empty):
      - NM comparison is primary determination.
      - Deletions (len(gref)>=1): genome fetch + ±10 bp refinement validates the Ref.
      - Insertions (gref=''): genome validation not applicable.
      - agreement_str: 'NM_validated'|'NM_mismatch'|'NM_corrected'|'NM_N/A'|'refalt_unresolved'

    winning_allele: 'A', 'B', or None (probe_only)
    winning_ts    : alignment dict or None (probe_only)
    Returns (ref_char, alt_char, agreement_str, final_pos).
    final_pos may be refined relative to the input for deletion markers.
    """
    allele_a = candidates_info["AlleleA"]
    allele_b = candidates_info["AlleleB"]
    is_indel = len(allele_a) != 1 or len(allele_b) != 1

    if is_indel:
        # NM comparison is the only determination method for indels
        nm_result = None
        if winning_allele is not None and winning_ts is not None:
            nm_result = determine_ref_alt(winning_allele, winning_ts,
                                           ts_aligns, candidates_info)
        if nm_result is None:
            return None, None, "refalt_unresolved", final_pos

        ref_char, alt_char = nm_result
        gref = ref_char  # longer or non-empty allele

        if gref == "":
            # Single-base insertion: NM assigned empty string as Ref, but the
            # genome base at final_pos may match alt_char, meaning NM got it
            # backwards (genome has the base → it is Ref, empty is Alt).
            if len(alt_char) == 1:
                try:
                    genome_base = fasta.fetch(chr_, final_pos - 1, final_pos).upper()
                    alt_fwd = strand_normalize(alt_char, strand)
                    if genome_base == alt_fwd:
                        # Genome confirms alt_char is actually Ref; swap.
                        return alt_char, ref_char, "NM_corrected", final_pos
                except (ValueError, KeyError):
                    pass
            return ref_char, alt_char, "NM_N/A", final_pos

        # Deletion: validate gref against genome, with ±10 bp refinement.
        # If the deletion sequence cannot be found, return "NM_mismatch" so the
        # marker passes through to qc_filter.py's design-conflict filter for removal.
        gref_fwd = reverse_complement(gref) if strand == "-" else gref
        refined = _refine_deletion_pos(fasta, chr_, final_pos, gref_fwd)
        if refined is not None:
            return ref_char, alt_char, "NM_validated", refined
        return ref_char, alt_char, "NM_mismatch", final_pos

    # ── SNP path ──────────────────────────────────────────────────────────────
    # Method 1: genome lookup (primary)
    genome_result = resolve_ref_from_genome(
        fasta, chr_, final_pos, allele_a, allele_b, strand
    )

    # Method 2: NM comparison (parallel; only when TopSeq aligned)
    nm_result = None
    nm_tied   = False
    if winning_allele is not None and winning_ts is not None:
        nr = determine_ref_alt(winning_allele, winning_ts, ts_aligns, candidates_info)
        if nr is None:
            nm_tied = True
        else:
            nm_result = nr

    if genome_result is not None:
        if winning_allele is None:
            # probe_only: no NM available
            return genome_result[0], genome_result[1], "NM_N/A", final_pos
        if nm_tied:
            return genome_result[0], genome_result[1], "NM_tied", final_pos
        if nm_result is None:
            return genome_result[0], genome_result[1], "NM_N/A", final_pos
        agreement = "NM_match" if genome_result == nm_result else "NM_unmatch"
        return genome_result[0], genome_result[1], agreement, final_pos

    # Genome lookup failed
    if nm_result is not None:
        return nm_result[0], nm_result[1], "NM_only", final_pos

    return None, None, "refalt_unresolved", final_pos


# ── DECISION COUNTERS ─────────────────────────────────────────────────────────

@dataclass
class DecisionCounters:
    """
    Accumulates per-marker decision counts throughout the pipeline and prints
    a structured summary table on completion.
    """
    total_loaded:            int = 0
    # alignment group counters (raw, before any filtering)
    align_gp1: int = 0
    align_gp2: int = 0
    align_gp3: int = 0
    align_gp4: int = 0
    align_gp5: int = 0
    align_unmapped: int = 0   # nothing aligned — exits before rescue logic
    # valid triple filtering
    valid_pair_found:        int = 0
    no_valid_pair:           int = 0
    # no-valid-triple: TopSeq rescue path
    final_topseq_only:            int = 0
    topseq_rescue_failed_softclip: int = 0   # SNP target in soft-clipped region
    topseq_rescue_refalt_unresolved: int = 0   # topseq_only winner with RefAlt=refalt_unresolved
    # no-valid-triple: probe rescue path
    final_probe_only:            int = 0
    final_probe_rescue_locus_unresolved: int = 0  # probe locus unresolved (tie=locus_unresolved)
    probe_rescue_refalt_unresolved:      int = 0  # probe Ref/Alt unresolved (tie=*, RefAlt=refalt_unresolved)
    probe_rescue_unmapped:               int = 0  # probe rescue failed → unmapped (defensive; should be 0 for gp5)
    topseq_locus_unresolved_no_probe_rescue: int = 0  # topseq locus unresolved → no probe rescue
    # unresolved sub-totals by coordinate role (for Final Output nesting)
    final_topseq_n_probe_unresolved: int = 0  # position_unresolved + refalt_unresolved_main
    final_topseq_only_unresolved:    int = 0  # topseq_locus_unresolved_no_probe_rescue + topseq_rescue_refalt_unresolved
    final_probe_only_unresolved:     int = 0  # final_probe_rescue_locus_unresolved + probe_rescue_refalt_unresolved
    # tie-resolution breakdown for successful rescues
    topseq_rescue_tie_unique:    int = 0
    topseq_rescue_tie_as:        int = 0
    topseq_rescue_tie_das:       int = 0
    topseq_rescue_tie_nm:        int = 0
    topseq_rescue_tie_scaffold:  int = 0
    probe_rescue_tie_unique:     int = 0
    probe_rescue_tie_as:         int = 0
    probe_rescue_tie_das:        int = 0
    probe_rescue_tie_nm:         int = 0
    probe_rescue_tie_scaffold:   int = 0
    # position resolution — topseq_n_probe tie breakdown (all 6 values)
    unique_position:              int = 0   # tie=unique
    topseq_n_probe_tie_as:        int = 0   # tie=AS_resolved
    topseq_n_probe_tie_das:       int = 0   # tie=dAS_resolved
    nm_position_resolved:         int = 0   # tie=NM_resolved
    topseq_n_probe_tie_coorddelta: int = 0  # tie=CoordDelta_resolved
    scaffold_resolved:             int = 0  # tie=scaffold_resolved
    position_unresolved:           int = 0  # tie=locus_unresolved
    # ref/alt determination
    ref_alt_ref_resolved:    int = 0
    refalt_unresolved_main:  int = 0   # tie=unique|*_resolved + RefAlt=refalt_unresolved (topseq_n_probe main path)
    ref_base_mismatch:       int = 0
    strand_agreement_unexpected: int = 0   # StrandAgreementAsExpected == False
    # MAPQ distribution (topseq_n_probe winners only)
    mapq_60:                 int = 0
    mapq_30_59:              int = 0
    mapq_1_29:               int = 0
    mapq_0:                  int = 0
    # CoordDelta distribution (topseq_n_probe winners only)
    coord_delta_0:           int = 0   # probe == cigar exactly
    coord_delta_1:           int = 0   # small discrepancy → probe used
    coord_delta_ge2:         int = 0   # large discrepancy (≥2 bp) → cigar used
    coord_delta_neg1:        int = 0   # CIGAR unavailable (SNP in soft clip)
    # CoordSource breakdown (topseq_n_probe winners only)
    coord_source_probe:      int = 0   # MapInfo = probe coord
    coord_source_cigar:      int = 0   # MapInfo = cigar coord (CoordDelta≥2 or indel)
    # final outcome counts
    final_mapped:            int = 0
    final_scaffold_resolved:      int = 0
    final_nm_position_resolved:   int = 0
    final_unmapped:               int = 0
    final_unresolved:             int = 0

    def format_summary(self, three_d=None) -> str:
        W = 60
        SEP = "\u2500" * W
        lines = []

        def row(label, n, indent=0):
            pad = "  " * indent
            lines.append(f"{pad}  {label:<{W - 2 - len(pad)}} {n:>8,}")

        def hdr(title):
            lines.append(f"\u2500\u2500 {title} {SEP[len(title) + 4:]}")

        # Direct-counter fallback for the topseq_n_probe tie breakdown.
        # Used when three_d is None (e.g. unit tests that call format_summary directly).
        # In production three_d is always provided and gives Chr≠0-only counts;
        # the fallback uses total-per-tie counters (includes any Chr=0 RefAlt=refalt_unresolved).
        _tnp_counter = {
            "unique":              self.unique_position,
            "AS_resolved":         self.topseq_n_probe_tie_as,
            "dAS_resolved":        self.topseq_n_probe_tie_das,
            "NM_resolved":         self.nm_position_resolved,
            "CoordDelta_resolved": self.topseq_n_probe_tie_coorddelta,
            "scaffold_resolved":   self.scaffold_resolved,
        }

        def _nonzero(anchor, tie):
            """Chr\u22600 count for (anchor, tie); falls back to per-tie counters if no 3D dict."""
            if three_d is not None:
                return three_d.get((anchor, tie), {}).get("NM_*", 0)
            if anchor == "topseq_n_probe":
                return _tnp_counter.get(tie, 0)
            return 0

        # \u2500\u2500 derived totals \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        # mapq_total == topseq_n_probe Chr\u22600 == valid_pair_found - position_unresolved - refalt_unresolved_main
        mapq_total        = self.mapq_60 + self.mapq_30_59 + self.mapq_1_29 + self.mapq_0
        topseq_only_total = self.final_topseq_only + self.final_topseq_only_unresolved
        probe_only_total  = self.final_probe_only  + self.final_probe_only_unresolved
        total_check       = (self.valid_pair_found + topseq_only_total
                             + probe_only_total + self.final_unmapped)

        lines.append("")
        lines.append("=== REMAPPING DECISION SUMMARY ===")
        lines.append(f"  {'Total markers loaded:':<{W - 2}} {self.total_loaded:>8,}")
        lines.append("")

        # \u2500\u2500 Step 1 \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        hdr("Step 1: Alignment Status (before filtering)")
        row("gp1 \u2013 both TopSeq alleles + probe:",   self.align_gp1)
        row("gp2 \u2013 one TopSeq allele  + probe:",    self.align_gp2)
        row("gp3 \u2013 both TopSeq alleles, no probe:", self.align_gp3)
        row("gp4 \u2013 one TopSeq allele, no probe:",   self.align_gp4)
        row("gp5 \u2013 probe only (no TopSeq):",        self.align_gp5)
        row("unmapped (nothing aligned):",               self.align_unmapped)
        lines.append("")

        # \u2500\u2500 Step 2: all three anchor branches as inline trees \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        hdr("Step 2: Valid Triple Filtering (chr + strand + overlap)")

        # \u2500\u2500 topseq_n_probe branch \u2500\u2500
        row("\u22651 valid triple \u2192 anchor=topseq_n_probe:", self.valid_pair_found)
        row("\u251c\u2500 tie=locus_unresolved (Chr=0):", self.position_unresolved, indent=1)
        row("\u251c\u2500 Ref/Alt=refalt_unresolved (Chr=0):", self.refalt_unresolved_main, indent=1)
        row("\u2514\u2500 Ref/Alt assigned \u2192 topseq_n_probe (Chr\u22600):", mapq_total, indent=1)
        for _t in ("unique", "AS_resolved", "dAS_resolved", "NM_resolved",
                   "CoordDelta_resolved", "scaffold_resolved"):
            row(f"\u2502  tie={_t}:", _nonzero("topseq_n_probe", _t), indent=3)

        # \u2500\u2500 no-valid-triple branches \u2500\u2500
        row("No valid triple (total):", self.no_valid_pair)
        lines.append("    TopSeq rescue (gp1\u2013gp4, anchor=topseq_only, coord=TopSeq CIGAR walk):")
        row("\u251c\u2500 tie=unique|*_resolved \u2192 Ref/Alt assigned (Chr\u22600):", self.final_topseq_only, indent=2)
        row("\u2502  tie=unique:",            self.topseq_rescue_tie_unique,   indent=3)
        row("\u2502  tie=AS_resolved:",       self.topseq_rescue_tie_as,       indent=3)
        row("\u2502  tie=dAS_resolved:",      self.topseq_rescue_tie_das,      indent=3)
        row("\u2502  tie=NM_resolved:",       self.topseq_rescue_tie_nm,       indent=3)
        row("\u2502  tie=scaffold_resolved:", self.topseq_rescue_tie_scaffold,  indent=3)
        row("\u251c\u2500 tie=unique|*_resolved + RefAlt=refalt_unresolved (Chr=0):", self.topseq_rescue_refalt_unresolved, indent=2)
        row("\u251c\u2500 tie=locus_unresolved (Chr=0):", self.topseq_locus_unresolved_no_probe_rescue, indent=2)
        row("\u2514\u2500 SNP in soft-clipped region \u2192 anchor=N/A (Chr=0):", self.topseq_rescue_failed_softclip, indent=2)
        lines.append("    Probe rescue (gp5, anchor=probe_only, coord=probe alignment):")
        row("\u251c\u2500 tie=unique|*_resolved \u2192 Ref/Alt assigned (Chr\u22600):", self.final_probe_only, indent=2)
        row("\u2502  tie=unique:",            self.probe_rescue_tie_unique,    indent=3)
        row("\u2502  tie=AS_resolved:",       self.probe_rescue_tie_as,        indent=3)
        row("\u2502  tie=dAS_resolved:",      self.probe_rescue_tie_das,       indent=3)
        row("\u2502  tie=NM_resolved:",       self.probe_rescue_tie_nm,        indent=3)
        row("\u2502  tie=scaffold_resolved:", self.probe_rescue_tie_scaffold,   indent=3)
        row("\u251c\u2500 tie=unique|*_resolved + RefAlt=refalt_unresolved (Chr=0):", self.probe_rescue_refalt_unresolved, indent=2)
        row("\u2514\u2500 tie=locus_unresolved (Chr=0):", self.final_probe_rescue_locus_unresolved, indent=2)
        lines.append("")

        # \u2500\u2500 Diagnostics \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        coord_diag_total = (self.coord_delta_0 + self.coord_delta_1
                            + self.coord_delta_ge2 + self.coord_delta_neg1)

        hdr(f"Diagnostics (topseq_n_probe Chr\u22600: {mapq_total:,})")
        lines.append(
            f"  = {self.valid_pair_found:,} total \u2212 {self.position_unresolved:,} "
            f"tie=locus_unresolved \u2212 {self.refalt_unresolved_main:,} Ref/Alt=refalt_unresolved"
        )
        lines.append("")
        lines.append("  MAPQ Distribution (min of winning TopSeq+probe pair):")
        row("MAPQ = 60:",   self.mapq_60,    indent=1)
        row("MAPQ 30\u201359:", self.mapq_30_59, indent=1)
        row("MAPQ  1\u201329:", self.mapq_1_29,  indent=1)
        row("MAPQ = 0:",    self.mapq_0,     indent=1)
        lines.append("")
        lines.append(f"  CoordDelta Distribution (of {coord_diag_total:,} markers):")
        row("CoordDelta = 0    (probe = cigar):",              self.coord_delta_0,    indent=1)
        row("CoordDelta = 1    (small diff \u2192 probe_cigar used):",  self.coord_delta_1,    indent=1)
        row("CoordDelta \u2265 2    (large diff \u2192 topseq_cigar used):", self.coord_delta_ge2,  indent=1)
        row("CoordDelta = \u22121   (SNP in soft clip, no topseq cigar):", self.coord_delta_neg1, indent=1)
        lines.append("")
        lines.append("  CoordSource Breakdown:")
        row("probe_cigar  (MapInfo = probe CIGAR coord):",  self.coord_source_probe, indent=1)
        row("topseq_cigar (MapInfo = TopSeq CIGAR coord):", self.coord_source_cigar, indent=1)
        lines.append("")
        lines.append("  Ref/Alt Diagnostics (topseq_n_probe Chr\u22600):")
        row("NM tie resolved by ref lookup:", self.ref_alt_ref_resolved, indent=1)
        row("Genome ref base mismatches (RefBaseMatch column):", self.ref_base_mismatch, indent=1)
        row("Strand agreement unexpected (StrandAgreementAsExpected=False):", self.strand_agreement_unexpected, indent=1)
        lines.append("")

        # \u2500\u2500 3D summary table \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        lines.append("\u2550" * W)
        lines.append("3-Dimension Summary (anchor \u00d7 tie \u00d7 Ref/Alt outcome)")
        if three_d is not None:
            ANCHOR_ORDER = ["topseq_n_probe", "topseq_only", "probe_only", "N/A"]
            TIE_ORDER    = ["unique", "AS_resolved", "dAS_resolved", "NM_resolved",
                            "CoordDelta_resolved", "scaffold_resolved", "locus_unresolved", "N/A"]

            lines.append(
                f"  {'anchor / tie':<28} {'NM_*(Chr\u22600)':>10}"
                f" {'unresolved(Chr=0)':>17} {'N/A(Chr=0)':>10} {'Total':>8}"
            )
            lines.append("  " + "\u2500" * 73)

            grand_nm = grand_amb = grand_na = 0
            for anchor in ANCHOR_ORDER:
                anchor_total = sum(
                    v
                    for (a, t), d in three_d.items() if a == anchor
                    for v in d.values()
                )
                if anchor_total == 0:
                    continue
                lines.append(f"  anchor={anchor}")
                for t in TIE_ORDER:
                    d  = three_d.get((anchor, t), {})
                    nm  = d.get("NM_*", 0)
                    amb = d.get("refalt_unresolved", 0)
                    na  = d.get("N/A", 0)
                    row_total = nm + amb + na
                    if row_total == 0:
                        continue
                    lines.append(
                        f"    tie={t:<24} {nm:>10,} {amb:>17,} {na:>10,} {row_total:>8,}"
                    )
                    grand_nm  += nm
                    grand_amb += amb
                    grand_na  += na
            grand_total = grand_nm + grand_amb + grand_na
            lines.append("  " + "\u2500" * 73)
            lines.append(
                f"  {'Total':<28} {grand_nm:>10,} {grand_amb:>17,}"
                f" {grand_na:>10,} {grand_total:>8,}"
            )
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
            # extract_candidates already normalises "-" → ""; no need to re-strip here.
            pre, a, b, post = extract_candidates(raw_topseq)
            if pre is not None:
                seq_a = pre + a + post
                seq_b = pre + b + post
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
    col_align_status = f"AlignmentStatus_{assembly}"
    col_anchor       = f"anchor_{assembly}"
    col_tie          = f"tie_{assembly}"
    col_refalt_agree = f"RefAltMethodAgreement_{assembly}"
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
        col_align_status: [],
        col_anchor:       [],
        col_tie:          [],
        col_refalt_agree: [],
    }
    # Side-file row lists (one per anchor). Each row records an unresolved marker's
    # competing-alignment detail. `UnresolvedReason` distinguishes locus-tie from
    # refalt-tie inside each file.
    tnp_unresolved_rows = []     # main path: topseq_n_probe
    to_unresolved_rows  = []     # rescue path: topseq_only
    po_unresolved_rows  = []     # rescue path: probe_only
    scaffold_rows  = []
    nm_pos_rows    = []
    counters       = DecisionCounters()
    counters.total_loaded = len(df)

    fasta = pysam.FastaFile(args.reference)
    try:

        def _append_unmapped_cols(anchor, tie="N/A"):
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
            new_cols[col_align_status].append("N/A")  # overwritten by caller when known
            new_cols[col_anchor].append(anchor)
            new_cols[col_tie].append(tie)
            new_cols[col_refalt_agree].append("N/A")

        for _, row in df.iterrows():
            name      = row["Name"]
            ts_aligns = raw_topseq.get(name, {})
            pb_aligns = probe_candidates.get(name, [])
            info      = candidates_info.get(name)

            # Alignment census (diagnostic — computed before any filtering)
            align_status = compute_alignment_status(ts_aligns, pb_aligns)

            # ΔScore: AS_best − AS_2nd across all TopSeq alignments.
            all_as = sorted(
                [a["AS"] for aligns in ts_aligns.values()
                 for a in aligns if a.get("AS", -1) >= 0],
                reverse=True,
            )
            delta_score = (all_as[0] - all_as[1]) if len(all_as) >= 2 else -1

            # Update alignment group counters
            if align_status == "gp1":   counters.align_gp1 += 1
            elif align_status == "gp2": counters.align_gp2 += 1
            elif align_status == "gp3": counters.align_gp3 += 1
            elif align_status == "gp4": counters.align_gp4 += 1
            elif align_status == "gp5": counters.align_gp5 += 1

            if not info or align_status == "unmapped":
                _append_unmapped_cols("N/A", "N/A")
                new_cols[col_align_status][-1] = align_status
                counters.align_unmapped += 1
                counters.final_unmapped += 1
                continue

            probe_seq   = str(row.get("AlleleA_ProbeSeq", ""))
            topseq_a    = info.get("TopSeqA", "")
            topseq_b    = info.get("TopSeqB", "")
            assay       = assay_types.get(name, "II")

            # ── Locus anchoring ──────────────────────────────────────────────
            triples = build_valid_triples(ts_aligns, pb_aligns,
                                           probe_seq, topseq_a, topseq_b)

            winning_allele = winning_ts = winning_pb = None
            anchor = tie_status = None

            if triples:
                counters.valid_pair_found += 1
                result = rank_and_resolve(triples, ts_aligns, pb_aligns,
                                           info, assay)
                if result[0] == "locus_unresolved":
                    _, competing = result
                    tnp_unresolved_rows.extend({"Name": name, **r} for r in competing)
                    _append_unmapped_cols("topseq_n_probe", "locus_unresolved")  # Case A
                    new_cols[col_align_status][-1] = align_status
                    counters.position_unresolved += 1
                    counters.final_unresolved += 1
                    counters.final_topseq_n_probe_unresolved += 1
                    continue

                tie_status     = result[0]
                winning_allele = result[1]
                winning_ts     = result[2]
                winning_pb     = result[3]
                anchor         = "topseq_n_probe"

                if tie_status == "scaffold_resolved":
                    scaffold_rows.extend({"Name": name, **r} for r in result[4])
                    counters.scaffold_resolved += 1
                elif tie_status in ("AS_resolved", "dAS_resolved",
                                    "NM_resolved", "CoordDelta_resolved"):
                    nm_pos_rows.extend({"Name": name, **r} for r in result[4])
                    if tie_status == "NM_resolved":
                        counters.nm_position_resolved += 1
                    elif tie_status == "AS_resolved":
                        counters.topseq_n_probe_tie_as += 1
                    elif tie_status == "dAS_resolved":
                        counters.topseq_n_probe_tie_das += 1
                    else:  # CoordDelta_resolved
                        counters.topseq_n_probe_tie_coorddelta += 1
                else:
                    counters.unique_position += 1

            else:
                # No valid triples — try TopSeq rescue
                counters.no_valid_pair += 1
                best_allele, best_ts, ts_tie, ts_pool = best_topseq_rescue(ts_aligns)

                if best_ts is not None and ts_tie != "locus_unresolved":
                    # TopSeq-only: derive coordinate from CIGAR
                    target_idx = info["PreLen"] if best_ts["Strand"] == "+" else info["PostLen"]
                    cigar_coord, cigar_in_sc = parse_cigar_to_ref_pos(
                        best_ts["Pos"], best_ts["Cigar"], target_idx
                    )
                    if cigar_in_sc or cigar_coord == 0:
                        _append_unmapped_cols("N/A", "N/A")
                        new_cols[col_align_status][-1] = align_status
                        counters.topseq_rescue_failed_softclip += 1
                        counters.final_unmapped += 1
                        continue

                    # Empty-allele CIGAR correction for deletions: when the winning
                    # allele's sequence is empty (deletion allele), seq[PreLen] lands
                    # on the first suffix base rather than the deletion start, placing
                    # cigar_coord len(other_allele) bases too far right.
                    is_indel = (info["AlleleA"] == "" or info["AlleleB"] == "")
                    if is_indel and not cigar_in_sc \
                            and info[f"Allele{best_allele}"] == "":
                        other_allele = "B" if best_allele == "A" else "A"
                        cigar_coord -= len(info[f"Allele{other_allele}"])

                    ref_alt_result = determine_ref_alt_v2(
                        best_allele, best_ts, ts_aligns, info,
                        fasta, best_ts["Chr"], cigar_coord, best_ts["Strand"]
                    )
                    if ref_alt_result[0] is None:
                        _append_unmapped_cols("topseq_only", ts_tie)  # Case D
                        new_cols[col_refalt_agree][-1] = "refalt_unresolved"
                        new_cols[col_align_status][-1] = align_status
                        counters.topseq_rescue_refalt_unresolved += 1
                        counters.final_unresolved += 1
                        counters.final_topseq_only_unresolved += 1
                        # Record the single winning alignment for the side file.
                        to_unresolved_rows.append({
                            "Name": name,
                            "UnresolvedReason": "topseq_refalt_tie",
                            "Rank": 1,
                            "TopSeqAllele": best_allele,
                            "Chr": best_ts["Chr"],
                            "Pos": best_ts["Pos"],
                            "Strand": best_ts["Strand"],
                            "MAPQ": best_ts["MAPQ"],
                            "NM": best_ts["NM"],
                            "AS": best_ts.get("AS", -1),
                        })
                        continue

                    ref_char, alt_char, refalt_agree, cigar_coord = ref_alt_result
                    ref_base_match_str = "N/A"
                    if len(ref_char) == 1 and len(alt_char) == 1:
                        try:
                            genome_base = fasta.fetch(
                                best_ts["Chr"], cigar_coord - 1, cigar_coord
                            ).upper()
                            ref_char_fwd = strand_normalize(ref_char, best_ts["Strand"])
                            ref_base_match_str = "True" if genome_base == ref_char_fwd else "False"
                            if ref_base_match_str == "False":
                                counters.ref_base_mismatch += 1
                        except (ValueError, KeyError):
                            ref_base_match_str = "False"
                            counters.ref_base_mismatch += 1

                    new_cols[col_chr].append(best_ts["Chr"])
                    new_cols[col_pos].append(cigar_coord)
                    new_cols[col_strand].append(best_ts["Strand"])
                    new_cols[col_ref].append(ref_char)
                    new_cols[col_alt].append(alt_char)
                    new_cols[col_mapq_ts].append(best_ts["MAPQ"])
                    new_cols[col_mapq_pb].append(float('nan'))
                    new_cols[col_delta].append(delta_score)
                    new_cols[col_qcov].append(compute_qcov(best_ts["Cigar"]))
                    new_cols[col_scfrac].append(compute_soft_clip_frac(best_ts["Cigar"]))
                    new_cols[col_cigar_coord].append(cigar_coord)
                    new_cols[col_probe_coord].append(0)
                    new_cols[col_coord_delta].append(-1)
                    new_cols[col_coord_source].append("topseq_cigar")
                    new_cols[col_ref_match].append(ref_base_match_str)
                    new_cols[col_probe_strand].append("N/A")
                    new_cols[col_strand_agree].append("N/A")
                    new_cols[col_align_status].append(align_status)
                    new_cols[col_anchor].append("topseq_only")
                    new_cols[col_tie].append(ts_tie)
                    new_cols[col_refalt_agree].append(refalt_agree)
                    counters.final_topseq_only += 1
                    if ts_tie == "unique":              counters.topseq_rescue_tie_unique   += 1
                    elif ts_tie == "AS_resolved":       counters.topseq_rescue_tie_as       += 1
                    elif ts_tie == "dAS_resolved":      counters.topseq_rescue_tie_das      += 1
                    elif ts_tie == "NM_resolved":       counters.topseq_rescue_tie_nm       += 1
                    elif ts_tie == "scaffold_resolved": counters.topseq_rescue_tie_scaffold += 1
                    continue

                else:
                    # TopSeq aligned but tie-break failed → probe cannot disambiguate.
                    # best_topseq_rescue returns (None, None, "locus_unresolved", pool)
                    # when TopSeq alignments were mapped but unresolvable — best_ts is
                    # None in that case, so check ts_tie as well.
                    if best_ts is not None or ts_tie == "locus_unresolved":
                        _append_unmapped_cols("topseq_only", "locus_unresolved")  # Case C
                        new_cols[col_align_status][-1] = align_status
                        counters.topseq_locus_unresolved_no_probe_rescue += 1
                        counters.final_unresolved += 1
                        counters.final_topseq_only_unresolved += 1
                        # Record competing-alignment detail from the locus-tie pool.
                        for rank, (allele_lbl, a) in enumerate(ts_pool or [], 1):
                            to_unresolved_rows.append({
                                "Name": name,
                                "UnresolvedReason": "topseq_locus_tie",
                                "Rank": rank,
                                "TopSeqAllele": allele_lbl,
                                "Chr": a["Chr"],
                                "Pos": a["Pos"],
                                "Strand": a["Strand"],
                                "MAPQ": a["MAPQ"],
                                "NM": a["NM"],
                                "AS": a.get("AS", -1),
                            })
                        continue

                    # TopSeq produced no alignments (gp5) — try probe-only rescue
                    best_pb, pb_tie, pb_pool = best_probe_rescue(pb_aligns)

                    if pb_tie == "locus_unresolved":
                        _append_unmapped_cols("probe_only", "locus_unresolved")  # Case E
                        new_cols[col_align_status][-1] = align_status
                        counters.final_probe_rescue_locus_unresolved += 1
                        counters.final_unresolved += 1
                        counters.final_probe_only_unresolved += 1
                        for rank, (_, a) in enumerate(pb_pool or [], 1):
                            po_unresolved_rows.append({
                                "Name": name,
                                "UnresolvedReason": "probe_locus_tie",
                                "Rank": rank,
                                "Chr": a["Chr"],
                                "Pos": a["Pos"],
                                "Strand": a["Strand"],
                                "MAPQ": a["MAPQ"],
                                "NM": a["NM"],
                                "AS": a.get("AS", -1),
                            })
                        continue

                    if best_pb is None:  # no mapped probe alignments (tie="N/A")
                        _append_unmapped_cols("N/A", pb_tie if pb_tie else "N/A")
                        new_cols[col_align_status][-1] = align_status
                        counters.probe_rescue_unmapped += 1
                        counters.final_unmapped += 1
                        continue

                    # Probe-only: derive coordinate from probe CIGAR
                    pb_coord = get_probe_coordinate(
                        best_pb["Pos"], best_pb["Cigar"], best_pb["Strand"], assay
                    )
                    ref_alt_result = determine_ref_alt_v2(
                        None, None, ts_aligns, info,
                        fasta, best_pb["Chr"], pb_coord, best_pb["Strand"]
                    )
                    if ref_alt_result[0] is None:
                        _append_unmapped_cols("probe_only", pb_tie)  # Case F
                        new_cols[col_refalt_agree][-1] = "refalt_unresolved"
                        new_cols[col_align_status][-1] = align_status
                        counters.probe_rescue_refalt_unresolved += 1
                        counters.final_unresolved += 1
                        counters.final_probe_only_unresolved += 1
                        po_unresolved_rows.append({
                            "Name": name,
                            "UnresolvedReason": "probe_refalt_tie",
                            "Rank": 1,
                            "Chr": best_pb["Chr"],
                            "Pos": best_pb["Pos"],
                            "Strand": best_pb["Strand"],
                            "MAPQ": best_pb["MAPQ"],
                            "NM": best_pb["NM"],
                            "AS": best_pb.get("AS", -1),
                        })
                        continue

                    ref_char, alt_char, refalt_agree, pb_coord = ref_alt_result

                    new_cols[col_chr].append(best_pb["Chr"])
                    new_cols[col_pos].append(pb_coord)
                    new_cols[col_strand].append(best_pb["Strand"])
                    new_cols[col_ref].append(ref_char)
                    new_cols[col_alt].append(alt_char)
                    new_cols[col_mapq_ts].append(float('nan'))
                    new_cols[col_mapq_pb].append(best_pb["MAPQ"])
                    new_cols[col_delta].append(-1)
                    new_cols[col_qcov].append(0.0)
                    new_cols[col_scfrac].append(0.0)
                    new_cols[col_cigar_coord].append(0)
                    new_cols[col_probe_coord].append(pb_coord)
                    new_cols[col_coord_delta].append(-1)
                    new_cols[col_coord_source].append("probe_cigar")
                    new_cols[col_ref_match].append("N/A")
                    new_cols[col_probe_strand].append(best_pb["Strand"])
                    new_cols[col_strand_agree].append("N/A")
                    new_cols[col_align_status].append(align_status)
                    new_cols[col_anchor].append("probe_only")
                    new_cols[col_tie].append(pb_tie)
                    new_cols[col_refalt_agree].append(refalt_agree)
                    counters.final_probe_only += 1
                    if pb_tie == "unique":              counters.probe_rescue_tie_unique   += 1
                    elif pb_tie == "AS_resolved":       counters.probe_rescue_tie_as       += 1
                    elif pb_tie == "dAS_resolved":      counters.probe_rescue_tie_das      += 1
                    elif pb_tie == "NM_resolved":       counters.probe_rescue_tie_nm       += 1
                    elif pb_tie == "scaffold_resolved": counters.probe_rescue_tie_scaffold += 1
                    continue

            # ── Winner path (topseq_n_probe) ─────────────────────────────────
            # Coordinate computation
            c_pos = get_probe_coordinate(
                winning_pb["Pos"], winning_pb["Cigar"],
                winning_pb["Strand"], assay
            )
            target_idx = info["PreLen"] if winning_ts["Strand"] == "+" else info["PostLen"]
            cigar_coord, cigar_in_sc = parse_cigar_to_ref_pos(
                winning_ts["Pos"], winning_ts["Cigar"], target_idx
            )
            is_indel = len(candidates_info[name]["AlleleA"]) != 1 or \
                       len(candidates_info[name]["AlleleB"]) != 1
            # For deletion markers where the winning allele is the empty (deletion)
            # allele, parse_cigar_to_ref_pos targets seq[PreLen] = SUFFIX[0], which
            # sits len(deleted_bases) past the true deletion-sequence start.
            # Subtract that length so cigar_coord points to the first deleted base.
            if is_indel and not cigar_in_sc \
                    and candidates_info[name][f"Allele{winning_allele}"] == "":
                other_allele = "B" if winning_allele == "A" else "A"
                cigar_coord -= len(candidates_info[name][f"Allele{other_allele}"])
            if cigar_in_sc:
                cigar_out       = 0
                coord_delta_val = -1
                final_pos       = c_pos
                coord_source    = "probe_cigar"
                counters.coord_delta_neg1 += 1
                counters.coord_source_probe += 1
            else:
                coord_delta_val = abs(c_pos - cigar_coord)
                cigar_out       = cigar_coord
                if coord_delta_val == 0:
                    counters.coord_delta_0 += 1
                elif coord_delta_val == 1:
                    counters.coord_delta_1 += 1
                else:
                    counters.coord_delta_ge2 += 1
                if is_indel or coord_delta_val >= 2:
                    final_pos    = cigar_coord
                    coord_source = "topseq_cigar"
                    counters.coord_source_cigar += 1
                else:
                    final_pos    = c_pos
                    coord_source = "probe_cigar"
                    counters.coord_source_probe += 1

            # Ref/Alt determination (uses final_pos = MapInfo)
            ref_alt_result = determine_ref_alt_v2(
                winning_allele, winning_ts, ts_aligns, info,
                fasta, winning_ts["Chr"], final_pos, winning_ts["Strand"]
            )
            if ref_alt_result[0] is None:
                tnp_unresolved_rows.append({
                    "Name": name,
                    "UnresolvedReason": "NM_and_genome_tie",
                    "PairRank": 1,
                    "TopSeqAllele": winning_allele,
                    "TopSeqChr": winning_ts["Chr"],
                    "TopSeqPos": winning_ts["Pos"],
                    "TopSeqStrand": winning_ts["Strand"],
                    "TopSeqMAPQ": winning_ts["MAPQ"],
                    "TopSeqNM": winning_ts["NM"],
                    "ProbeChr": winning_pb["Chr"],
                    "ProbePos": winning_pb["Pos"],
                    "ProbeMAPQ": winning_pb["MAPQ"],
                    "MinMAPQ": min(winning_ts["MAPQ"], winning_pb["MAPQ"]),
                })
                _append_unmapped_cols("topseq_n_probe", tie_status)  # Case B
                new_cols[col_refalt_agree][-1] = "refalt_unresolved"
                new_cols[col_align_status][-1] = align_status
                counters.refalt_unresolved_main += 1
                counters.final_unresolved += 1
                counters.final_topseq_n_probe_unresolved += 1
                continue

            ref_char, alt_char, refalt_agree, final_pos = ref_alt_result

            # Deletion minus-strand coordinate correction
            if len(ref_char) > len(alt_char) and winning_ts["Strand"] == "-":
                if coord_source == "probe_cigar":
                    final_pos -= len(ref_char) - len(alt_char)
                    c_pos     -= len(ref_char) - len(alt_char)

            # Probe strand agreement (sequence-derived; reporting only)
            pb_strand, sa_expected = compute_probe_strand_agreement(
                topseq_strand=winning_ts["Strand"],
                probe_align_strand=winning_pb["Strand"],
                probe_seq=probe_seq,
                topseq_a=topseq_a,
                topseq_b=topseq_b,
            )
            if sa_expected == "False":
                counters.strand_agreement_unexpected += 1

            # RefBaseMatch
            ref_base_match_str = "N/A"
            if len(ref_char) == 1 and len(alt_char) == 1:
                try:
                    genome_base = fasta.fetch(
                        winning_ts["Chr"], final_pos - 1, final_pos
                    ).upper()
                    ref_char_fwd = strand_normalize(ref_char, winning_ts["Strand"])
                    ref_base_match_str = "True" if genome_base == ref_char_fwd else "False"
                    if ref_base_match_str == "False":
                        counters.ref_base_mismatch += 1
                except (ValueError, KeyError):
                    ref_base_match_str = "False"
                    counters.ref_base_mismatch += 1

            # MAPQ distribution counter (kept for report)
            min_mapq = min(winning_ts["MAPQ"], winning_pb["MAPQ"])
            if min_mapq == 60:   counters.mapq_60    += 1
            elif min_mapq >= 30: counters.mapq_30_59 += 1
            elif min_mapq >= 1:  counters.mapq_1_29  += 1
            else:                counters.mapq_0     += 1

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
            new_cols[col_probe_coord].append(c_pos)
            new_cols[col_coord_delta].append(coord_delta_val)
            new_cols[col_coord_source].append(coord_source)
            new_cols[col_ref_match].append(ref_base_match_str)
            new_cols[col_probe_strand].append(pb_strand)
            new_cols[col_strand_agree].append(sa_expected)
            new_cols[col_align_status].append(align_status)
            new_cols[col_anchor].append(anchor)
            new_cols[col_tie].append(tie_status)
            new_cols[col_refalt_agree].append(refalt_agree)

            if tie_status == "scaffold_resolved":
                counters.final_scaffold_resolved += 1
            elif tie_status == "NM_resolved":
                counters.final_nm_position_resolved += 1
            elif refalt_agree in ("NM_only", "NM_tied"):
                counters.ref_alt_ref_resolved += 1
                counters.final_mapped += 1
            else:
                counters.final_mapped += 1

        # ── Step 6: Write outputs ────────────────────────────────────────────────
        for col, vals in new_cols.items():
            df[col] = vals

        print(f"[remap] Writing output: {args.output}")
        df.to_csv(args.output, index=False)

        # Three side-output files, one per anchor. Each file pools all unresolved
        # markers for that path (locus-tie + RefAlt-tie); the UnresolvedReason
        # column distinguishes them inside the file.
        for label, rows in (
            ("topseq_n_probe_unresolved_markers.csv", tnp_unresolved_rows),
            ("topseq_only_unresolved_markers.csv",    to_unresolved_rows),
            ("probe_only_unresolved_markers.csv",     po_unresolved_rows),
        ):
            if rows:
                p = os.path.join(out_dir, label)
                pd.DataFrame(rows).to_csv(p, index=False)
                print(f"[remap] Unresolved markers written to: {p}")

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
        # Compute 3D groupby (anchor × tie × RefAlt bucket) from final columns
        three_d = defaultdict(lambda: defaultdict(int))
        for i in range(len(new_cols[col_anchor])):
            a  = new_cols[col_anchor][i]
            t  = new_cols[col_tie][i]
            r  = str(new_cols[col_refalt_agree][i])
            # Buckets: "NM_*" (resolved) / "refalt_unresolved" (Chr=0) / "N/A"
            rb = "NM_*" if r.startswith("NM_") else r
            three_d[(a, t)][rb] += 1

        summary = counters.format_summary(three_d)
        print(summary)
        report_path = os.path.join(out_dir, "remapping_Report.txt")
        with open(report_path, "w") as f:
            f.write(summary + "\n")
        print(f"[remap] Remapping report: {report_path}")

    finally:
        fasta.close()


if __name__ == "__main__":
    run_remapping(parse_args())
