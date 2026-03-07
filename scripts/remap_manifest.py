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
  MAPQ_Probe              Mapping quality of the selected probe alignment (0 = fallback used)

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

import pandas as pd

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
    return p.parse_args()

# ── MANIFEST PARSING ─────────────────────────────────────────────────────────

def locate_data_section(filename, start_marker="[Assay]", end_marker="[Controls]"):
    """Returns (skiprows, nrows) for the data section of an Illumina manifest."""
    header_line = 0
    footer_line = None
    with open(filename) as f:
        for i, line in enumerate(f):
            stripped = line.strip()
            if stripped == start_marker:
                header_line = i + 1
            if stripped == end_marker:
                footer_line = i
                break
    nrows = (footer_line - header_line) if footer_line is not None else None
    return header_line, nrows

# ── SEQUENCE HELPERS ─────────────────────────────────────────────────────────

def extract_candidates(top_seq):
    """
    Parses TopGenomicSeq format 'PREFIX[A/B]SUFFIX' into (pre, alleleA, alleleB, post).
    Returns (None, None, None, None) if the format is not recognised.
    """
    m = re.search(r"(.*?)\[(.*?)/(.*?)\](.*)", top_seq or "")
    if m:
        return m.group(1), m.group(2), m.group(3), m.group(4)
    return None, None, None, None

# ── CIGAR UTILITIES ──────────────────────────────────────────────────────────

def _parse_cigar(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def cigar_ref_span(cigar):
    """Reference bases consumed by the alignment (for computing alignment end position)."""
    return sum(n for n, op in _parse_cigar(cigar) if op in "MDN=X")


def get_alignment_end(pos, cigar):
    """1-based end coordinate of an alignment (inclusive)."""
    return pos + cigar_ref_span(cigar) - 1


def parse_cigar_to_ref_pos(start_pos, cigar, query_index):
    """
    Maps a 0-based query index to a 1-based reference coordinate via CIGAR.
    Used as the TopGenomicSeq fallback when no probe overlap is found.
    """
    ops = _parse_cigar(cigar)
    curr_q = 0
    curr_r = start_pos
    for n, op in ops:
        if op in "M=X":
            if query_index < curr_q + n:
                return curr_r + (query_index - curr_q)
            curr_q += n
            curr_r += n
        elif op in "IS":
            if query_index < curr_q + n:
                return curr_r  # inside insertion/clip: return junction
            curr_q += n
        elif op in "DN":
            curr_r += n
    return curr_r


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

# ── SAM PARSING ──────────────────────────────────────────────────────────────

def _get_nm(cols):
    for tag in cols[11:]:
        if tag.startswith("NM:i:"):
            return int(tag.split(":")[2])
    return 999


def parse_topseq_sam(sam_path):
    """
    Reads the TopGenomicSeq SAM and returns a dict:
      { snp_name: { 'A': {NM, Chr, Pos, Cigar, MAPQ, Strand, End}, 'B': ... } }
    Only primary, mapped alignments are kept.
    """
    results = {}
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@") or line.startswith("[M"):
                continue
            cols = line.split("\t")
            flag = int(cols[1])
            if flag & 4:   # unmapped
                continue
            if flag & 256: # secondary
                continue
            qname_full = cols[0]
            which = qname_full[-1]          # 'A' or 'B'
            qname = qname_full[:-2]         # strip '_A' or '_B'
            pos = int(cols[3])
            cigar = cols[5]
            entry = {
                "NM": _get_nm(cols),
                "Chr": cols[2],
                "Pos": pos,
                "Cigar": cigar,
                "MAPQ": int(cols[4]),
                "Strand": "-" if flag & 16 else "+",
                "End": get_alignment_end(pos, cigar),
            }
            results.setdefault(qname, {})[which] = entry
    return results


def parse_probe_sam(sam_path):
    """
    Reads the probe SAM and returns a dict:
      { snp_name: [ {Chr, Pos, Cigar, Strand, MAPQ, End}, ... ] }
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
                "End": get_alignment_end(pos, cigar),
            }
            results.setdefault(cols[0], []).append(entry)
    return results

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

            pre, a, b, post = extract_candidates(row.get("TopGenomicSeq", ""))
            if pre is not None:
                seq_a = pre + (a if a != "-" else "") + post
                seq_b = pre + (b if b != "-" else "") + post
                ft.write(f">{name}_A\n{seq_a}\n")
                ft.write(f">{name}_B\n{seq_b}\n")
                candidates_info[name] = {
                    "PreLen": len(pre), "PostLen": len(post),
                    "AlleleA": a, "AlleleB": b,
                }

    # ── Align with minimap2 ──────────────────────────────────────────────────
    print("[remap] Aligning sequences with minimap2...")
    t = args.threads
    subprocess.check_call(
        f"minimap2 -ax sr -t {t} {args.reference} {topseq_fasta} > {topseq_sam} 2>/dev/null",
        shell=True,
    )
    subprocess.check_call(
        # -N 5: allow up to 5 secondary alignments for overlap checking
        f"minimap2 -ax sr -N 5 -t {t} {args.reference} {probe_fasta} > {probe_sam} 2>/dev/null",
        shell=True,
    )

    # ── Parse TopGenomicSeq alignments → Ref/Alt authority ──────────────────
    print("[remap] Parsing TopGenomicSeq alignments...")
    raw_top = parse_topseq_sam(topseq_sam)
    topseq_results = {}

    for name, res in raw_top.items():
        if name not in candidates_info:
            continue
        info = candidates_info[name]
        res_a = res.get("A")
        res_b = res.get("B")
        nm_a = res_a["NM"] if res_a else 999
        nm_b = res_b["NM"] if res_b else 999

        # Lower edit distance = closer to reference → that allele is the Ref
        if nm_b < nm_a:
            winner, final_ref, final_alt = res_b, info["AlleleB"], info["AlleleA"]
        elif res_a:
            winner, final_ref, final_alt = res_a, info["AlleleA"], info["AlleleB"]
        else:
            winner, final_ref, final_alt = res_b, info["AlleleB"], info["AlleleA"]

        if winner:
            topseq_results[name] = {
                "Ref": final_ref, "Alt": final_alt,
                "Chr": winner["Chr"],
                "Start": winner["Pos"], "End": winner["End"],
                "MAPQ": winner["MAPQ"],
                "Cigar": winner["Cigar"],
                "Strand": winner["Strand"],
                "PreLen": info["PreLen"], "PostLen": info["PostLen"],
            }

    # ── Parse probe alignments ───────────────────────────────────────────────
    print("[remap] Parsing probe alignments...")
    probe_candidates = parse_probe_sam(probe_sam)

    # ── Resolve final coordinates ────────────────────────────────────────────
    print("[remap] Resolving coordinates...")
    new_cols = {
        col_chr: [], col_pos: [], col_strand: [],
        col_ref: [], col_alt: [],
        col_mapq_ts: [], col_mapq_pb: [],
    }

    for _, row in df.iterrows():
        name = row["Name"]
        ts = topseq_results.get(name)

        if not ts:
            # Failed to align TopGenomicSeq at all
            new_cols[col_chr].append("0")
            new_cols[col_pos].append(0)
            new_cols[col_strand].append("N/A")
            new_cols[col_ref].append("N")
            new_cols[col_alt].append("N")
            new_cols[col_mapq_ts].append(0)
            new_cols[col_mapq_pb].append(0)
            continue

        c_ref, c_alt   = ts["Ref"], ts["Alt"]
        c_chr, c_strand = ts["Chr"], ts["Strand"]
        c_mapq_ts = ts["MAPQ"]
        c_mapq_pb = 0

        # Select the probe alignment with the greatest overlap with the TopSeq region,
        # using MAPQ as a tiebreaker when overlaps are equal.
        best_pb = None
        best_overlap = 0
        best_mapq = -1
        for pb in probe_candidates.get(name, []):
            if pb["Chr"] != c_chr:
                continue
            ov = calculate_overlap(ts["Start"], ts["End"], pb["Pos"], pb["End"])
            if ov > best_overlap or (ov == best_overlap and pb["MAPQ"] > best_mapq):
                best_overlap = ov
                best_mapq = pb["MAPQ"]
                best_pb = pb

        if best_pb and best_pb["Strand"] == c_strand:
            # Case A: probe is consistent with TopSeq strand → high-precision coordinate
            c_mapq_pb = best_pb["MAPQ"]
            assay = assay_types.get(name, "II")
            raw_pos = get_probe_coordinate(best_pb["Pos"], best_pb["Cigar"], c_strand, assay)

            # Deletion correction on minus strand:
            # The probe's 3' end points to the high coordinate of the deletion event;
            # subtract the deletion length to get the canonical VCF position.
            if len(c_ref) > len(c_alt) and c_strand == "-":
                raw_pos -= (len(c_ref) - len(c_alt))

            c_pos = raw_pos
        else:
            # Case B: no consistent probe → fall back to TopGenomicSeq CIGAR
            target_idx = ts["PreLen"] if c_strand == "+" else ts["PostLen"]
            c_pos = parse_cigar_to_ref_pos(ts["Start"], ts["Cigar"], target_idx)

        new_cols[col_chr].append(c_chr)
        new_cols[col_pos].append(c_pos)
        new_cols[col_strand].append(c_strand)
        new_cols[col_ref].append(c_ref)
        new_cols[col_alt].append(c_alt)
        new_cols[col_mapq_ts].append(c_mapq_ts)
        new_cols[col_mapq_pb].append(c_mapq_pb)

    # ── Append columns and write output ─────────────────────────────────────
    for col, vals in new_cols.items():
        df[col] = vals

    print(f"[remap] Writing output: {args.output}")
    df.to_csv(args.output, index=False)

    print(f"[remap] Temp files written to: {temp_dir}")
    print("[remap] Done.")


if __name__ == "__main__":
    run_remapping(parse_args())
