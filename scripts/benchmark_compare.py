#!/usr/bin/env python3
"""
benchmark_compare.py — Compare remap_manifest.py output against a ground-truth
EquCab3-native Illumina manifest. Produces per-marker TSVs and a summary report.

Usage:
    python scripts/benchmark_compare.py \
        --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
        --remapped  results_E80selv2_to_equCab3/..._remapped_equCab3.csv \
        --assembly  equCab3 \
        [--output-dir results_E80selv2_to_equCab3/benchmark/] \
        [--baseline  results_E80selv2_to_equCab3/benchmark/benchmark_TIMESTAMP.tsv]
"""

import argparse
import os
import re
import sys
from datetime import datetime

import pandas as pd

from _strand_utils import strand_normalize, complement

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Benchmark remap_manifest.py output against ground-truth manifest.")
    p.add_argument("--manifest",    required=True, help="Illumina manifest CSV (EquCab3-native, GenomeBuild=3)")
    p.add_argument("--remapped",    required=True, help="*_remapped_{assembly}.csv from remap_manifest.py")
    p.add_argument("--assembly",    default=None,
                   help="Assembly label (e.g. equCab3). Optional — auto-detected from "
                        "the remapped CSV header if omitted. Pass explicitly only when "
                        "the CSV contains multiple assemblies.")
    p.add_argument("--reference",   default=None,
                   help="Reference genome FASTA (enables the explanatory layer with context-based verdicts)")
    p.add_argument("--output-dir",  default="./benchmark_out", help="Output directory (default: ./benchmark_out)")
    p.add_argument("--baseline",    default=None,  help="Path to prior benchmark_<timestamp>.tsv for diff")
    p.add_argument("--flank-len",   type=int, default=20,
                   help="Flank length for context-match checks (default: 20 bp)")
    p.add_argument("--traced",      default=None,
                   help="Path to {prefix}_remapped_{assembly}_traced.csv from qc_filter.py; "
                        "enables the QC filtration impact section")
    return p.parse_args()


# ── ASSEMBLY AUTO-DETECTION ───────────────────────────────────────────────────

_ANCHOR_SUFFIX_RE = re.compile(r"^anchor_(.+)$")


def detect_assembly(remapped_path: str, explicit=None) -> str:
    """Return the assembly label embedded in the remapped CSV's column names.

    If *explicit* is truthy, returns it unchanged (user override).
    Otherwise scans the CSV header for `anchor_<X>` columns and:
      - exactly one match → returns X
      - zero matches     → raises ValueError (not a current-schema CSV)
      - multiple matches → raises ValueError (ambiguous; ask user to pass --assembly)
    """
    if explicit:
        return explicit
    header = pd.read_csv(remapped_path, nrows=0)
    suffixes = set()
    for col in header.columns:
        m = _ANCHOR_SUFFIX_RE.match(col)
        if m:
            suffixes.add(m.group(1))
    if len(suffixes) == 1:
        return suffixes.pop()
    if not suffixes:
        raise ValueError(
            f"Remapped CSV {remapped_path!r} has no anchor_<assembly> column. "
            "Re-run remap_manifest.py to produce the current-schema output."
        )
    raise ValueError(
        f"Remapped CSV {remapped_path!r} contains multiple assemblies: "
        f"{sorted(suffixes)!r}. Ambiguous; pass --assembly explicitly."
    )


# ── CHROMOSOME NORMALISATION ──────────────────────────────────────────────────

def normalise_chr(chrom: str) -> str:
    """Normalise X chromosome aliases to 'X'. All other values pass through unchanged."""
    if str(chrom).startswith("X_"):
        return "X"
    return str(chrom)


# ── ALLELE PARSING ────────────────────────────────────────────────────────────

_BRACKET_RE = re.compile(r"\[(.+?)/(.+?)\]")
AMBIGUOUS_PAIRS = (frozenset({"A", "T"}), frozenset({"C", "G"}))


def parse_snp_alleles(snp_col):
    """Parse manifest SNP column (e.g. '[A/G]') into (A, B) single-base pair.

    Returns None for indels ('[D/I]'), empty/None input, or malformed strings.
    """
    if not snp_col:
        return None
    m = _BRACKET_RE.search(str(snp_col))
    if not m:
        return None
    a, b = m.group(1), m.group(2)
    # Reject indel-label markers
    if {a.upper(), b.upper()} == {"D", "I"}:
        return None
    # Require single-base alleles; multi-base means indel-in-SNP-column, not a SNP
    if len(a) != 1 or len(b) != 1:
        return None
    return a.upper(), b.upper()


def parse_topseq_alleles(topseq_col):
    """Parse TopGenomicSeq 'PREFIX[A/B]SUFFIX' into (prefix, A, B, suffix).

    The '-' placeholder (as in '[CCC/-]') is normalised to empty string so that
    deletion / insertion alleles are represented by their real sequences.
    Returns None for malformed input.
    """
    if not topseq_col:
        return None
    m = _BRACKET_RE.search(str(topseq_col))
    if not m:
        return None
    pre_end = m.start()
    post_start = m.end()
    a = "" if m.group(1) == "-" else m.group(1)
    b = "" if m.group(2) == "-" else m.group(2)
    return str(topseq_col)[:pre_end], a, b, str(topseq_col)[post_start:]


def alleles_match_snp(remapped_ref, remapped_alt, remapped_strand,
                      manifest_pair) -> bool:
    """True iff the remapped SNP alleles match the manifest {A,B} pair.

    Both strand-normalised and complemented sets are tried to accommodate SNP
    columns that use the opposite strand convention. Returns False for indels
    (manifest_pair is None).
    """
    if manifest_pair is None:
        return False
    if not remapped_ref or not remapped_alt:
        return False
    fwd_ref = strand_normalize(remapped_ref, remapped_strand).upper()
    fwd_alt = strand_normalize(remapped_alt, remapped_strand).upper()
    fwd_set = {fwd_ref, fwd_alt}
    manifest_set = {manifest_pair[0].upper(), manifest_pair[1].upper()}
    if fwd_set == manifest_set:
        return True
    comp_set = {complement(fwd_ref).upper(), complement(fwd_alt).upper()}
    return comp_set == manifest_set


# ── CONTEXT CHECKING (explanatory layer) ──────────────────────────────────────

def check_flanking_context(fasta, chrom, mapinfo, prefix, suffix,
                            allele_len, flank_len=20):
    """Return (forward_match, reverse_match) for the genome flanking at mapinfo.

    Forward test: genome[left] == prefix[-flank:] AND genome[right] == suffix[:flank].
    Reverse test: genome[left] == RC(suffix)[-flank:] AND genome[right] == RC(prefix)[:flank].

    mapinfo is 1-based start of the variant base; allele_len is the length of the
    reference allele at the locus (0 for a pure insertion; L for a deletion / SNP).
    """
    if not prefix or not suffix:
        return (False, False)
    k = min(flank_len, len(prefix), len(suffix))
    if k == 0:
        return (False, False)
    left_start = mapinfo - 1 - k
    if left_start < 0:
        return (False, False)
    right_start = mapinfo - 1 + allele_len
    try:
        left  = fasta.fetch(chrom, left_start, mapinfo - 1).upper()
        right = fasta.fetch(chrom, right_start, right_start + k).upper()
    except (ValueError, KeyError):
        return (False, False)
    fwd = (left == prefix[-k:].upper()) and (right == suffix[:k].upper())
    rc_prefix = strand_normalize(prefix, "-").upper()
    rc_suffix = strand_normalize(suffix, "-").upper()
    rev = (left == rc_suffix[-k:]) and (right == rc_prefix[:k])
    return (fwd, rev)


# ── EXPLANATORY VERDICT ───────────────────────────────────────────────────────

def classify_explanatory(row) -> str:
    """Produce a single verdict string from raw explanatory signals.

    Expected keys in *row*: context_forward, context_reverse, manifest_strand,
    remapped_strand, coord_ok, is_ambiguous_snp, is_probe_only, is_indel,
    deletion_seq_ok, insertion_absent, result.

    The "no context match either strand" branch splits three ways based on
    the top-level `result` classification so the verdict names describe
    what actually happened:
      - result == "unmapped"   → pipeline_unmapped   (no position assigned)
      - result == "wrong_chr"  → pipeline_wrong_chr  (placed on wrong chromosome)
      - otherwise              → pipeline_wrong_locus (placed on right chromosome,
                                                       but position's context fails)
    """
    fwd         = bool(row.get("context_forward"))
    rev         = bool(row.get("context_reverse"))
    coord_ok    = bool(row.get("coord_ok"))
    m_strand    = row.get("manifest_strand")
    r_strand    = row.get("remapped_strand")
    is_amb      = bool(row.get("is_ambiguous_snp"))
    is_po       = bool(row.get("is_probe_only"))
    result      = row.get("result", "")

    if is_amb and fwd:
        return "ambiguous_snp"
    if fwd and (r_strand != m_strand):
        return "manifest_strand_wrong"
    if rev and (r_strand != m_strand) and not fwd:
        return "pipeline_wrong_strand"
    if fwd and not coord_ok:
        return "manifest_coord_wrong"
    if not fwd and not rev:
        if is_po:
            return "probe_only_inconclusive"
        if result == "unmapped":
            return "pipeline_unmapped"
        if result == "wrong_chr":
            return "pipeline_wrong_chr"
        return "pipeline_wrong_locus"
    return "unresolved"


# ── COMPARISON ────────────────────────────────────────────────────────────────

def classify_marker(manifest_row: dict, remapped_row: dict) -> str:
    """
    Classify a single marker into one of 6 outcome categories.

    Category priority (highest first):
      unmapped         — remapped Chr=0 or Strand=N/A
      locus_unresolved — remapped_status=="locus_unresolved" (pipeline rejected all candidate loci)
      wrong_chr        — Chr mismatch
      coord_off        — Chr match but MapInfo differs
      coord_correct_strand_wrong — Chr+MapInfo match but Strand differs
      correct          — all three match

    Strand comparison: when `remapped_probe_strand` is provided it takes priority
    over `remapped_strand` (TopSeq alignment strand is uncorrelated with
    RefStrand; probe alignment strand is the right ground-truth comparison).
    A probe strand of "N/A" (topseq_only markers without a probe alignment)
    exempts the marker from the strand check. Missing probe strand falls back
    to `remapped_strand` for backward compatibility.
    """
    m_chr    = normalise_chr(manifest_row["manifest_chr"])
    m_pos    = manifest_row["manifest_pos"]
    m_strand = manifest_row["manifest_strand"]

    r_chr    = normalise_chr(remapped_row["remapped_chr"])
    r_pos    = remapped_row["remapped_pos"]
    r_strand = remapped_row["remapped_strand"]
    r_probe_strand = remapped_row.get("remapped_probe_strand")
    r_status = remapped_row["remapped_status"]

    # Priority 1: unmapped
    if r_chr == "0" or r_strand == "N/A":
        return "unmapped"

    # Priority 2: locus_unresolved (note: this short-circuits all coordinate
    # checks, so a marker here with a wrong chromosome is still counted as
    # 'locus_unresolved', not 'wrong_chr')
    if r_status == "locus_unresolved":
        return "locus_unresolved"

    # Priority 3: wrong chromosome
    if r_chr != m_chr:
        return "wrong_chr"

    # Priority 4: coordinate off
    try:
        if int(r_pos) != int(m_pos):
            return "coord_off"
    except (ValueError, TypeError):
        return "coord_off"

    # Priority 5: strand wrong (probe strand preferred; N/A exempts; None falls back)
    if r_probe_strand in ("+", "-"):
        strand_to_compare = r_probe_strand
    elif r_probe_strand in ("N/A", "nan"):
        strand_to_compare = None          # topseq_only: exempt from strand check
    else:
        strand_to_compare = r_strand      # missing column: fall back to TopSeq strand

    if strand_to_compare is not None and str(strand_to_compare) != str(m_strand):
        return "coord_correct_strand_wrong"

    return "correct"


# ── MANIFEST LOADING ──────────────────────────────────────────────────────────

def _locate_assay_section(path: str) -> tuple[int, int | None]:
    """
    Returns (header_line, footer_line) as 0-based line indices.
    header_line  : the line containing the column names (line after [Assay])
    footer_line  : the line containing [Controls], or None if absent
    """
    header_line = None
    footer_line = None
    with open(path) as f:
        for i, line in enumerate(f):
            stripped = line.strip().rstrip(",")
            if stripped == "[Assay]":
                header_line = i + 1
            if stripped == "[Controls]":
                footer_line = i
                break
    if header_line is None:
        raise ValueError(f"[Assay] section not found in {path}")
    return header_line, footer_line


def load_manifest(path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load an Illumina manifest CSV and partition markers into three DataFrames:
      main_df  — autosome or X (benchmarked)
      chry_df  — Chr == Y
      chr0_df  — Chr == 0 or MapInfo == 0

    All DataFrames have columns:
      Name, manifest_chr, manifest_pos, manifest_strand, SNP, TopGenomicSeq
    Chr is normalised (X_NC_009175.3 → X) in all three.

    manifest_strand is read directly from RefStrand (± relative to reference).
    SourceStrand / IlmnStrand are NOT consulted — TOP/BOT/PLUS/MINUS do not map to
    reference-strand orientation reliably.
    """
    header_line, footer_line = _locate_assay_section(path)
    nrows = (footer_line - header_line - 1) if footer_line is not None else None

    df = pd.read_csv(
        path,
        skiprows=header_line,
        nrows=nrows,
        dtype={"Chr": str, "MapInfo": str},
        usecols=["Name", "Chr", "MapInfo", "RefStrand", "SNP", "TopGenomicSeq"],
        low_memory=False,
    )
    df = df.rename(columns={"Chr": "manifest_chr", "MapInfo": "manifest_pos",
                             "RefStrand": "manifest_strand"})
    df["manifest_chr"] = df["manifest_chr"].apply(normalise_chr)
    df["manifest_pos"] = pd.to_numeric(df["manifest_pos"], errors="coerce").fillna(0).astype(int)
    # RefStrand is already '+' or '-'; pass through unchanged.
    df["manifest_strand"] = df["manifest_strand"].astype(str)

    # Partition
    is_y    = df["manifest_chr"] == "Y"
    is_chr0 = (df["manifest_chr"] == "0") | (df["manifest_pos"] == 0)
    is_main = ~is_y & ~is_chr0

    return df[is_main].reset_index(drop=True), \
           df[is_y].reset_index(drop=True), \
           df[is_chr0].reset_index(drop=True)


# ── REMAPPED LOADING ──────────────────────────────────────────────────────────

def load_remapped(path: str, assembly: str) -> pd.DataFrame:
    """
    Load the remapped CSV produced by remap_manifest.py.
    Returns a DataFrame with columns:
      Name, remapped_chr, remapped_pos, remapped_strand, remapped_status
    and optionally coord_delta (if CoordDelta_{assembly} is present in the CSV).
    Chr is normalised (X aliases → X).

    Requires the current-schema columns `anchor_{assembly}` and `tie_{assembly}`
    (produced by remap_manifest.py's 3-dimension output framework). Legacy
    CSVs that only carry the pre-3-D `MappingStatus_{assembly}` column are
    not supported — re-run remap_manifest.py to upgrade.

    remapped_status derivation:
      anchor == "topseq_only"              → "topseq_only"
      tie    == "locus_unresolved"         → "locus_unresolved"
      anchor == "N/A" (unmapped)           → "unmapped"
      otherwise                            → "mapped"
    """
    col_chr    = f"Chr_{assembly}"
    col_pos    = f"MapInfo_{assembly}"
    col_strand = f"Strand_{assembly}"
    col_probe_strand = f"ProbeStrand_{assembly}"
    col_ref    = f"Ref_{assembly}"
    col_alt    = f"Alt_{assembly}"
    col_anchor = f"anchor_{assembly}"
    col_tie    = f"tie_{assembly}"
    col_delta  = f"CoordDelta_{assembly}"

    header = pd.read_csv(path, nrows=0)

    required = ["Name", col_chr, col_pos, col_strand, col_anchor, col_tie]
    missing = [c for c in required if c not in header.columns]
    if missing:
        raise ValueError(
            f"Remapped CSV {path!r} is missing expected columns: {missing}\n"
            f"Re-run remap_manifest.py against this manifest to produce the "
            f"current-schema output (anchor_{assembly} + tie_{assembly})."
        )
    has_delta = col_delta in header.columns
    has_refalt = col_ref in header.columns and col_alt in header.columns
    has_probe_strand = col_probe_strand in header.columns

    usecols = ["Name", col_chr, col_pos, col_strand, col_anchor, col_tie]
    if has_delta:
        usecols.append(col_delta)
    if has_refalt:
        usecols += [col_ref, col_alt]
    if has_probe_strand:
        usecols.append(col_probe_strand)

    df = pd.read_csv(
        path,
        dtype={col_chr: str},
        usecols=usecols,
        low_memory=False,
    )
    if df["Name"].duplicated().any():
        dups = df["Name"][df["Name"].duplicated()].head(5).tolist()
        raise ValueError(
            f"Remapped CSV {path!r} contains duplicate Name values "
            f"(first few: {dups}). Each marker must have a unique Name."
        )

    # Build a unified remapped_status column
    def _derive_status(row):
        anchor = row[col_anchor]
        tie    = row[col_tie]
        # pandas reads "N/A" as NaN; treat NaN anchor as unmapped
        anchor_str = "" if pd.isna(anchor) else str(anchor)
        tie_str    = "" if pd.isna(tie)    else str(tie)
        if anchor_str == "topseq_only":
            return "topseq_only"
        if tie_str == "locus_unresolved":
            return "locus_unresolved"
        if anchor_str in ("N/A", ""):
            return "unmapped"
        return "mapped"
    df["remapped_status"] = df.apply(_derive_status, axis=1)
    df = df.rename(columns={col_anchor: "anchor", col_tie: "tie"})

    rename_map = {
        col_chr:    "remapped_chr",
        col_pos:    "remapped_pos",
        col_strand: "remapped_strand",
    }
    if has_delta:
        rename_map[col_delta] = "coord_delta"
    if has_refalt:
        rename_map[col_ref] = "remapped_ref"
        rename_map[col_alt] = "remapped_alt"
    if has_probe_strand:
        rename_map[col_probe_strand] = "remapped_probe_strand"
    df = df.rename(columns=rename_map)
    df["remapped_chr"] = df["remapped_chr"].apply(normalise_chr)
    df["remapped_pos"] = pd.to_numeric(df["remapped_pos"], errors="coerce")
    return df


# ── COMPARISON ────────────────────────────────────────────────────────────────

def compare_all(manifest_df: pd.DataFrame, remapped_df: pd.DataFrame,
                 fasta=None, flank_len: int = 20) -> pd.DataFrame:
    """
    Left-join manifest onto remapped, classify each marker, compute coord_offset.
    Markers present in manifest but absent from remapped are classified as 'unmapped'.

    When *fasta* (pysam.FastaFile) is supplied and the manifest includes
    TopGenomicSeq/SNP columns, the output also contains per-marker explanatory
    signals (context_forward, context_reverse, deletion_seq_ok, insertion_absent,
    is_ambiguous_snp, is_probe_only, allele_ok) and a synthesised verdict column.

    Returns a DataFrame with all manifest columns plus:
      remapped_chr, remapped_pos, remapped_strand, remapped_status,
      result, coord_offset,
      (optionally) the explanatory columns listed above + verdict.
    """
    merged = manifest_df.merge(remapped_df, on="Name", how="left")

    # Fill missing remapped rows (not in remapped CSV at all)
    merged["remapped_chr"]    = merged["remapped_chr"].fillna("0")
    merged["remapped_pos"]    = merged["remapped_pos"].fillna(0)
    merged["remapped_strand"] = merged["remapped_strand"].fillna("N/A")
    merged["remapped_status"] = merged["remapped_status"].fillna("unmapped")
    if "remapped_probe_strand" in merged.columns:
        merged["remapped_probe_strand"] = merged["remapped_probe_strand"].fillna("N/A")
    if "anchor" in merged.columns:
        merged["anchor"] = merged["anchor"].fillna("N/A")
    if "tie" in merged.columns:
        merged["tie"] = merged["tie"].fillna("N/A")

    results = []
    offsets = []
    for _, row in merged.iterrows():
        m = {k: row[k] for k in ("manifest_chr", "manifest_pos", "manifest_strand")}
        r = {k: row[k] for k in ("remapped_chr", "remapped_pos",
                                  "remapped_strand", "remapped_status")}
        if "remapped_probe_strand" in merged.columns:
            r["remapped_probe_strand"] = row["remapped_probe_strand"]
        cat = classify_marker(m, r)
        results.append(cat)
        if cat == "coord_off":
            try:
                offsets.append(int(row["remapped_pos"]) - int(row["manifest_pos"]))
            except (ValueError, TypeError):
                offsets.append(None)
        else:
            offsets.append(None)

    merged["result"]       = results
    merged["coord_offset"] = offsets
    if "coord_delta" in merged.columns:
        merged["coord_delta"] = merged["coord_delta"].fillna(-1).astype(int)

    # Explanatory-layer signals (only when fasta + manifest context are available)
    if fasta is not None and "TopGenomicSeq" in merged.columns:
        _add_explanatory_signals(merged, fasta, flank_len)

    return merged


def _add_explanatory_signals(merged: pd.DataFrame, fasta, flank_len: int) -> None:
    """Populate context_forward, context_reverse, allele_ok, deletion_seq_ok,
    insertion_absent, is_ambiguous_snp, is_probe_only, verdict columns in place.

    Only runs detailed signal computation on non-`correct` rows; passing markers
    receive empty / default values.
    """
    cf_list, cr_list = [], []
    ao_list, dso_list, ia_list = [], [], []
    amb_list, po_list = [], []
    verdicts = []

    has_refalt = "remapped_ref" in merged.columns and "remapped_alt" in merged.columns
    # Vectorised probe-only mask avoids a per-row label→positional lookup over 80k+ rows.
    if "anchor" in merged.columns:
        probe_only_mask = (merged["anchor"] == "probe_only").to_numpy()
    else:
        probe_only_mask = [False] * len(merged)

    for pos_i, (idx, row) in enumerate(merged.iterrows()):
        result = row["result"]
        is_probe_only = bool(probe_only_mask[pos_i])
        po_list.append(is_probe_only)

        snp_pair = parse_snp_alleles(row.get("SNP"))
        is_amb = snp_pair is not None and frozenset(snp_pair) in AMBIGUOUS_PAIRS
        amb_list.append(is_amb)

        ts_parsed = parse_topseq_alleles(row.get("TopGenomicSeq"))

        # allele_ok for SNPs (when we have remapped_ref/remapped_alt)
        allele_ok = None
        if snp_pair is not None and has_refalt:
            r_ref = row.get("remapped_ref")
            r_alt = row.get("remapped_alt")
            if pd.notna(r_ref) and pd.notna(r_alt) and row["remapped_strand"] in ("+", "-"):
                allele_ok = alleles_match_snp(str(r_ref), str(r_alt),
                                               row["remapped_strand"], snp_pair)
        ao_list.append(allele_ok)

        # Default signals (used only for non-correct rows)
        cf = cr = False
        dso = ia = None

        if result != "correct" and ts_parsed is not None and row["remapped_strand"] in ("+", "-"):
            pre, allele_a, allele_b, post = ts_parsed
            # Determine allele length at the locus.
            # For SNPs both alleles are length 1. For indels one allele is empty.
            # Use max(len(a), len(b)) to skip over any deletion sequence when checking the
            # right-side flank, and 0 when it's an insertion (ref has no allele bases).
            if allele_a == "" or allele_b == "":
                # deletion: ref allele is the non-empty one; insertion: ref allele is empty
                ref_is_deletion = (allele_b == "" and allele_a != "") or (allele_a != "" and len(allele_a) > len(allele_b))
                allele_len = len(allele_a if allele_a != "" else allele_b) if ref_is_deletion else 0
            else:
                allele_len = 1

            try:
                mapinfo = int(row["remapped_pos"])
                if mapinfo > 0 and row["remapped_chr"] not in ("0", "", None):
                    cf, cr = check_flanking_context(
                        fasta, row["remapped_chr"], mapinfo, pre, post,
                        allele_len, flank_len=flank_len,
                    )
                    # Indel-specific signals
                    deletion_allele = allele_a if (allele_a and not allele_b) else (allele_b if (allele_b and not allele_a) else None)
                    if deletion_allele is not None and len(deletion_allele) > 0:
                        try:
                            fetched = fasta.fetch(row["remapped_chr"], mapinfo - 1,
                                                   mapinfo - 1 + len(deletion_allele)).upper()
                            fwd_ok = (fetched == deletion_allele.upper())
                            rev_ok = (fetched == strand_normalize(deletion_allele, "-").upper())
                            dso = (fwd_ok or rev_ok)
                        except (ValueError, KeyError):
                            dso = False
                    insertion_allele = None
                    if allele_a == "" and allele_b != "":
                        insertion_allele = allele_b
                    elif allele_b == "" and allele_a != "":
                        insertion_allele = allele_a if (allele_a != "" and len(allele_a) < len(allele_b or "")) else None
                    if insertion_allele is not None and len(insertion_allele) > 0:
                        try:
                            fetched = fasta.fetch(row["remapped_chr"], mapinfo - 1,
                                                   mapinfo - 1 + len(insertion_allele)).upper()
                            ia = (fetched != insertion_allele.upper())
                        except (ValueError, KeyError):
                            ia = None
            except (ValueError, TypeError):
                pass

        cf_list.append(cf)
        cr_list.append(cr)
        dso_list.append(dso)
        ia_list.append(ia)

        # Compute verdict (empty string for `correct` rows). Probe strand is the
        # authoritative strand for the RefStrand comparison.
        #   "+"/"-"  → use as r_strand_for_verdict
        #   "N/A"    → topseq_only; exempt (use manifest strand so "disagree" branch is false)
        #   missing  → fall back to TopSeq strand
        if result == "correct":
            verdicts.append("")
        else:
            probe_strand = row.get("remapped_probe_strand")
            if probe_strand in ("+", "-"):
                r_strand_for_verdict = probe_strand
            elif probe_strand == "N/A":
                r_strand_for_verdict = row["manifest_strand"]   # force-equal so strand branches skip
            else:
                r_strand_for_verdict = row["remapped_strand"]
            v_row = {
                "context_forward": cf,
                "context_reverse": cr,
                "manifest_strand": row["manifest_strand"],
                "remapped_strand": r_strand_for_verdict,
                "coord_ok": (result not in ("coord_off", "wrong_chr", "unmapped")),
                "is_ambiguous_snp": is_amb,
                "is_probe_only": is_probe_only,
                "is_indel": (ts_parsed is not None and (ts_parsed[1] == "" or ts_parsed[2] == "")),
                "deletion_seq_ok": dso,
                "insertion_absent": ia,
                "result": result,
            }
            verdicts.append(classify_explanatory(v_row))

    merged["context_forward"]  = cf_list
    merged["context_reverse"]  = cr_list
    merged["allele_ok"]        = ao_list
    merged["deletion_seq_ok"]  = dso_list
    merged["insertion_absent"] = ia_list
    merged["is_ambiguous_snp"] = amb_list
    merged["is_probe_only"]    = po_list
    merged["verdict"]          = verdicts


# ── OUTPUT ────────────────────────────────────────────────────────────────────

_TSV_COLS = [
    "Name",
    "manifest_chr", "manifest_pos", "manifest_strand",
    "remapped_chr",  "remapped_pos",  "remapped_strand", "remapped_probe_strand",
    "remapped_status", "anchor", "tie",
    "result", "coord_offset", "coord_delta",
    # Explanatory layer (present when --reference was supplied)
    "allele_ok",
    "context_forward", "context_reverse",
    "deletion_seq_ok", "insertion_absent",
    "is_ambiguous_snp", "is_probe_only",
    "verdict",
    # QC filtration trace (present when --traced was supplied)
    "why_filtered",
]


# QC filter stage labels — must match scripts/qc_filter.py:WHY_FILTERED_LABELS.
# Redefined here so benchmark_compare.py doesn't import from qc_filter.
QC_STAGE_ORDER = [
    "stage_1_failed_markers",
    "stage_2_design_conflict",
    "stage_3_min_anchor",
    "stage_4_tie_policy",
    "stage_5_min_refalt_confidence",
    "stage_6_mapq_topseq",
    "stage_7_mapq_probe",
    "stage_8_coord_delta",
    "stage_9_indel_excluded",
    "stage_10_polymorphic",
    "stage_11_ambiguous_snp",
]

CATEGORIES = [
    "correct",
    "coord_correct_strand_wrong",
    "coord_off",
    "wrong_chr",
    "unmapped",
    "locus_unresolved",
]


def _compute_transitions(
    curr: pd.Series, base: pd.Series
) -> tuple[int, dict[tuple[str, str], int]]:
    """Return (changed_count, transitions_dict) between two result series indexed by Name."""
    common = curr.index.intersection(base.index)
    changed = common[curr[common] != base[common]]
    transitions: dict[tuple[str, str], int] = {}
    for name in changed:
        key = (base[name], curr[name])
        transitions[key] = transitions.get(key, 0) + 1
    return len(changed), transitions


def write_tsv(df: pd.DataFrame, path: str) -> None:
    """Write result DataFrame to a tab-separated file using the standard column order."""
    cols = [c for c in _TSV_COLS if c in df.columns]
    df[cols].to_csv(path, sep="\t", index=False)


def load_baseline(path: str) -> pd.DataFrame:
    """Load a previously written benchmark TSV for diff comparison."""
    return pd.read_csv(path, sep="\t", dtype={"manifest_chr": str, "remapped_chr": str})


def _fmt(n: int, total: int) -> str:
    pct = 100 * n / total if total else 0.0
    return f"{n:>10,}  ({pct:5.1f}%)"


def _marker_type(name: str) -> str:
    if str(name).startswith("Affx-"):
        return "AFFX"
    if "ilmndup" in str(name):
        return "ilmndup"
    return "standard"


def stratify_by_coord_delta(result_df: pd.DataFrame) -> list[dict] | None:
    """
    Stratify accuracy by CoordDelta bucket.

    CoordDelta = |probe_coord - CIGAR_coord|; -1 means CIGAR coord was unavailable
    (SNP target fell in a soft-clipped region).

    Returns a list of dicts (one per non-empty bucket), each with:
      label, n, coord_accurate (correct + coord_correct_strand_wrong), correct

    Returns None if the coord_delta column is absent (old remapped CSVs).
    """
    if "coord_delta" not in result_df.columns:
        return None

    delta = result_df["coord_delta"]
    buckets = [
        ("delta = 0",    delta == 0),
        ("delta = 1",    delta == 1),
        ("delta = 2-10", (delta >= 2) & (delta <= 10)),
        ("delta > 10",   delta > 10),
        ("delta = -1",   delta == -1),
    ]

    rows = []
    for label, mask in buckets:
        subset = result_df[mask]
        n = len(subset)
        if n == 0:
            continue
        n_correct       = int((subset["result"] == "correct").sum())
        n_strand_wrong  = int((subset["result"] == "coord_correct_strand_wrong").sum())
        rows.append({
            "label":          label,
            "n":              n,
            "coord_accurate": n_correct + n_strand_wrong,
            "correct":        n_correct,
        })
    return rows


def build_three_d_accuracy(result_df: pd.DataFrame) -> dict | None:
    """
    Group result_df by (anchor, tie) and count each benchmark result category.

    Returns a dict mapping (anchor_str, tie_str) → {category: count, ...},
    or None if anchor/tie columns are absent (old remapped CSVs).
    """
    if "anchor" not in result_df.columns or "tie" not in result_df.columns:
        return None
    three_d = {}
    for (anchor, tie), group in result_df.fillna("N/A").groupby(["anchor", "tie"]):
        three_d[(str(anchor), str(tie))] = {
            cat: int((group["result"] == cat).sum()) for cat in CATEGORIES
        }
    return three_d


def format_three_d_accuracy_table(three_d: dict) -> str:
    """
    Format the 3-dimension accuracy breakdown as a string matching the style of
    qc_filter.format_three_d_table.

    Columns: correct | coord_off | wrong_chr | other | Total | Acc%
      other = unmapped + locus_unresolved + coord_correct_strand_wrong
      Acc%  = (correct + coord_correct_strand_wrong) / Total  (coord-accurate)
    """
    ANCHOR_ORDER = ["topseq_n_probe", "topseq_only", "probe_only", "N/A"]
    TIE_ORDER    = ["unique", "AS_resolved", "dAS_resolved", "NM_resolved",
                    "CoordDelta_resolved", "scaffold_resolved", "locus_unresolved", "N/A"]
    W = 88

    lines = [
        "═" * W,
        "3-Dimension Accuracy Breakdown  (anchor × tie × benchmark outcome)",
        f"  {'anchor / tie':<28} {'correct':>10} {'coord_off':>10}"
        f" {'wrong_chr':>10} {'other':>7} {'Total':>8} {'Acc%':>7}",
        "  " + "─" * 84,
    ]

    grand = {cat: 0 for cat in CATEGORIES}

    for anchor in ANCHOR_ORDER:
        anchor_data = {t: d for (a, t), d in three_d.items() if a == anchor}
        if not anchor_data or sum(v for d in anchor_data.values() for v in d.values()) == 0:
            continue
        lines.append(f"  anchor={anchor}")
        for tie in TIE_ORDER:
            d = anchor_data.get(tie)
            if d is None:
                continue
            total = sum(d.values())
            if total == 0:
                continue
            correct    = d["correct"]
            strand_bad = d["coord_correct_strand_wrong"]
            coord_off  = d["coord_off"]
            wrong_chr  = d["wrong_chr"]
            other      = d["unmapped"] + d["locus_unresolved"] + strand_bad
            acc        = 100.0 * (correct + strand_bad) / total
            lines.append(
                f"    tie={tie:<24} {correct:>10,} {coord_off:>10,}"
                f" {wrong_chr:>10,} {other:>7,} {total:>8,} {acc:>6.1f}%"
            )
            for cat in CATEGORIES:
                grand[cat] += d[cat]

    grand_total  = sum(grand.values())
    grand_strand = grand["coord_correct_strand_wrong"]
    grand_acc    = 100.0 * (grand["correct"] + grand_strand) / grand_total if grand_total else 0.0
    grand_other  = grand["unmapped"] + grand["locus_unresolved"] + grand_strand
    lines += [
        "  " + "─" * 84,
        f"  {'Total':<28} {grand['correct']:>10,} {grand['coord_off']:>10,}"
        f" {grand['wrong_chr']:>10,} {grand_other:>7,} {grand_total:>8,} {grand_acc:>6.1f}%",
    ]
    return "\n".join(lines)


# ── QC FILTRATION IMPACT ──────────────────────────────────────────────────────

def load_traced_why_filtered(path: str, assembly: str) -> pd.DataFrame:
    """Load `Name` and `WhyFiltered_{assembly}` from the traced CSV.

    Returns a DataFrame with columns: Name, why_filtered.
    Markers that passed all filters have why_filtered == "" (pandas reads empty
    CSV fields as NaN; we coerce those back to "").
    """
    col = f"WhyFiltered_{assembly}"
    df = pd.read_csv(path, usecols=["Name", col], low_memory=False, dtype=str)
    df = df.rename(columns={col: "why_filtered"})
    df["why_filtered"] = df["why_filtered"].fillna("")
    return df


def compute_qc_impact(result_df: pd.DataFrame, stages=QC_STAGE_ORDER) -> dict:
    """Compute QC-impact metrics. *result_df* must have 'result' and 'why_filtered'.

    Keys returned:
      passed_n, passed_correct, passed_accuracy_pct
      confusion       : DataFrame indexed by why_filtered, columns = benchmark categories
      per_stage       : list of (stage, n_removed, n_non_correct, precision_pct)
      cumulative      : list of (after_stage_label, n_remaining, n_correct, accuracy_pct)
      fp_df           : markers with result='correct' but why_filtered != ''  (removed by QC but actually correct)
      fn_df           : markers with why_filtered == '' but result != 'correct' (kept by QC but non-correct)
    """
    passed = result_df[result_df["why_filtered"] == ""]
    passed_n = len(passed)
    passed_correct = int((passed["result"] == "correct").sum())
    passed_acc = 100.0 * passed_correct / passed_n if passed_n else 0.0

    # Confusion matrix — rows: filter-stage label (or "" for passed); cols: benchmark category
    confusion = pd.crosstab(result_df["why_filtered"], result_df["result"])

    # Per-stage precision — fraction of markers removed by the stage that were non-correct
    per_stage = []
    for stage in stages:
        sub = result_df[result_df["why_filtered"] == stage]
        n = len(sub)
        if n == 0:
            continue
        nc = int((sub["result"] != "correct").sum())
        per_stage.append((stage, n, nc, 100.0 * nc / n))

    # Cumulative passing-set accuracy — starting from all benchmarked, remove each stage in order
    running = result_df.copy()
    cumulative = [(
        "(before QC)",
        len(running),
        int((running["result"] == "correct").sum()),
        100.0 * (running["result"] == "correct").sum() / len(running) if len(running) else 0.0,
    )]
    for stage in stages:
        running = running[running["why_filtered"] != stage]
        n = len(running)
        c = int((running["result"] == "correct").sum())
        acc = 100.0 * c / n if n else 0.0
        cumulative.append((stage, n, c, acc))

    fp_df = result_df[(result_df["why_filtered"] != "") & (result_df["result"] == "correct")]
    fn_df = result_df[(result_df["why_filtered"] == "") & (result_df["result"] != "correct")]

    return {
        "passed_n":              passed_n,
        "passed_correct":        passed_correct,
        "passed_accuracy_pct":   passed_acc,
        "confusion":             confusion,
        "per_stage":             per_stage,
        "cumulative":            cumulative,
        "fp_df":                 fp_df,
        "fn_df":                 fn_df,
    }


def format_qc_impact_section(impact: dict) -> list[str]:
    """Render `compute_qc_impact` output as a list of report lines."""
    lines = []
    W = 68

    def w(s=""):
        lines.append(s)

    w("QC FILTRATION IMPACT")
    w(f"  Markers surviving all QC filters: {impact['passed_n']:>8,}")
    w(f"    of which correct:               {impact['passed_correct']:>8,}  ({impact['passed_accuracy_pct']:5.1f}%)")
    w()

    # Confusion matrix
    w("  Confusion matrix  (stage × benchmark result)")
    conf = impact["confusion"]
    cats = [c for c in ["correct", "coord_correct_strand_wrong", "coord_off",
                         "wrong_chr", "unmapped", "locus_unresolved"] if c in conf.columns]
    header = f"    {'stage':<28} " + " ".join(f"{c[:10]:>10}" for c in cats) + f" {'total':>8}"
    w(header)
    w("    " + "-" * (len(header) - 4))
    # "" (passed) row first
    if "" in conf.index:
        row = conf.loc[""]
        row_total = int(row.sum())
        row_str = " ".join(f"{int(row.get(c, 0)):>10,}" for c in cats)
        w(f"    {'(passed all QC)':<28} {row_str} {row_total:>8,}")
    # then each stage that has removals
    for stage in QC_STAGE_ORDER:
        if stage not in conf.index:
            continue
        row = conf.loc[stage]
        row_total = int(row.sum())
        row_str = " ".join(f"{int(row.get(c, 0)):>10,}" for c in cats)
        w(f"    {stage:<28} {row_str} {row_total:>8,}")
    w()

    # Per-stage precision
    w("  Per-stage precision  (% of removed that are non-correct)")
    w(f"    {'stage':<28} {'removed':>10} {'non-correct':>12} {'precision':>10}")
    w("    " + "-" * 64)
    for stage, n, nc, pct in impact["per_stage"]:
        w(f"    {stage:<28} {n:>10,} {nc:>12,} {pct:>9.1f}%")
    w()

    # Cumulative accuracy
    w("  Cumulative passing-set accuracy  (stages applied in order)")
    w(f"    {'after':<28} {'remaining':>10} {'correct':>10} {'accuracy':>10}")
    w("    " + "-" * 62)
    for label, n, c, acc in impact["cumulative"]:
        w(f"    {label:<28} {n:>10,} {c:>10,} {acc:>9.1f}%")
    w()

    # False positives: correct markers removed by QC
    fp = impact["fp_df"]
    w(f"  False positives  (correct markers removed by QC): {len(fp):,}")
    if len(fp) > 0:
        fp_by_stage = fp["why_filtered"].value_counts()
        w("    FP by stage:")
        for stage, n in fp_by_stage.items():
            w(f"      {stage:<28} {int(n):>6,}")
        w()
        w(f"    First {min(10, len(fp))} FP markers:")
        w(f"      {'Name':<35} {'stage':<28} {'Chr':>4} {'MapInfo':>12}")
        for _, row in fp.head(10).iterrows():
            w(f"      {str(row.get('Name',''))[:35]:<35} "
              f"{str(row.get('why_filtered',''))[:28]:<28} "
              f"{str(row.get('manifest_chr','')):>4} "
              f"{str(row.get('manifest_pos','')):>12}")
        w()

    # False negatives: non-correct markers surviving QC
    fn = impact["fn_df"]
    w(f"  False negatives  (non-correct markers surviving QC): {len(fn):,}")
    if len(fn) > 0:
        w("    FN by benchmark result:")
        for cat, n in fn["result"].value_counts().items():
            w(f"      {cat:<28} {int(n):>6,}")
        if "verdict" in fn.columns:
            w("    FN by explanatory verdict:")
            vc = fn["verdict"].fillna("").value_counts()
            for verdict, n in vc.items():
                label = verdict if verdict else "(no verdict)"
                w(f"      {label:<28} {int(n):>6,}")
        w()

    return lines


def write_report(
    result_df: pd.DataFrame,
    chry_df: pd.DataFrame,
    chr0_df: pd.DataFrame,
    path: str,
    assembly: str,
    manifest_path: str,
    remapped_path: str,
    baseline_df: pd.DataFrame | None = None,
) -> None:
    """Write the human-readable benchmark report."""
    total_manifest = len(result_df) + len(chry_df) + len(chr0_df)
    benchmarked    = len(result_df)
    lines = []

    def w(s=""):
        lines.append(s)

    w(f"Benchmark Report — assembly: {assembly}")
    w(f"Run: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    w(f"Manifest:  {manifest_path}")
    w(f"Remapped:  {remapped_path}")
    w("-" * 60)
    w()
    w("SCOPE")
    w(f"  Total manifest markers:        {total_manifest:>10,}")
    w(f"  Benchmarked (autosome + X):    {benchmarked:>10,}")
    w(f"  Excluded Chr=Y:                {len(chry_df):>10,}")
    w(f"  Excluded Chr=0:                {len(chr0_df):>10,}")
    w()
    w("HEADLINE COUNTS")
    for cat in CATEGORIES:
        n = (result_df["result"] == cat).sum()
        w(f"  {cat:<32} {_fmt(n, benchmarked)}")
    w()
    w("BREAKDOWN BY MARKER TYPE")
    for mtype in ("standard", "AFFX", "ilmndup"):
        subset = result_df[result_df["Name"].apply(_marker_type) == mtype]
        if len(subset) == 0:
            continue
        w(f"  {mtype} ({len(subset):,} markers):")
        for cat in CATEGORIES:
            n = (subset["result"] == cat).sum()
            if n > 0:
                w(f"    {cat:<30} {_fmt(n, len(subset))}")
    w()

    # Coordinate offset distribution (coord_off markers only)
    coord_off = result_df[result_df["result"] == "coord_off"].copy()
    if len(coord_off) > 0:
        w("COORDINATE OFFSET DISTRIBUTION  (coord_off markers only)")
        offsets = coord_off["coord_offset"].dropna().abs()
        buckets = [
            ("offset = 1 bp",                      offsets == 1),
            ("offset = 2–10 bp",                   (offsets >= 2) & (offsets <= 10)),
            ("offset = 11–50 bp",                  (offsets >= 11) & (offsets <= 50)),
            ("offset = 51 bp (probe-strand bug)",  offsets == 51),
            ("offset = 52+ bp",                    offsets >= 52),
        ]
        for label, mask in buckets:
            w(f"  {label:<40} {mask.sum():>8,}")
        w()

    # CoordDelta stratification (only when column present)
    delta_rows = stratify_by_coord_delta(result_df)
    if delta_rows is not None:
        w("ACCURACY STRATIFIED BY COORD_DELTA")
        w("  CoordDelta = |probe_coord - CIGAR_coord|; -1 = CIGAR target in soft clip")
        w("  coord-accurate = correct + coord_correct_strand_wrong")
        w()
        w(f"  {'bucket':<14}  {'N':>8}  {'coord-accurate':>20}  {'correct':>20}")
        w("  " + "-" * 68)
        for r in delta_rows:
            ca_str = _fmt(r["coord_accurate"], r["n"])
            co_str = _fmt(r["correct"],        r["n"])
            w(f"  {r['label']:<14}  {r['n']:>8,}  {ca_str}  {co_str}")
        w()

    # 3D accuracy breakdown (only when anchor/tie columns are present)
    three_d = build_three_d_accuracy(result_df)
    if three_d is not None:
        w("3-DIMENSION ACCURACY BREAKDOWN")
        w("  Acc% = (correct + coord_correct_strand_wrong) / Total  [coord-accurate]")
        w("  other = unmapped + locus_unresolved + coord_correct_strand_wrong")
        w()
        w(format_three_d_accuracy_table(three_d))
        w()

    # Verdict distribution (explanatory layer; present when --reference was used)
    if "verdict" in result_df.columns:
        non_correct = result_df[result_df["result"] != "correct"]
        if len(non_correct) > 0:
            w("EXPLANATORY VERDICTS  (non-correct markers only)")
            counts = non_correct["verdict"].value_counts()
            for verdict, n in counts.items():
                w(f"  {str(verdict):<30} {_fmt(int(n), len(non_correct))}")
            w()

    # QC filtration impact (present when --traced was supplied)
    if "why_filtered" in result_df.columns:
        w("-" * 60)
        impact = compute_qc_impact(result_df)
        for line in format_qc_impact_section(impact):
            w(line)

    # Diff section
    if baseline_df is not None:
        w("-" * 60)
        w("DIFF VS BASELINE")
        curr = result_df.set_index("Name")["result"]
        base = baseline_df.set_index("Name")["result"]
        n_changed, transitions = _compute_transitions(curr, base)
        w(f"  Markers that changed category: {n_changed:>8,}")
        w()
        for (from_cat, to_cat), count in sorted(transitions.items(), key=lambda x: -x[1]):
            w(f"  {from_cat:<32} → {to_cat:<32} {count:>8,}")
        w()

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_diff(result_df: pd.DataFrame, baseline_df: pd.DataFrame, path: str) -> None:
    """Write a standalone diff file listing category transitions vs baseline."""
    curr = result_df.set_index("Name")["result"]
    base = baseline_df.set_index("Name")["result"]
    n_changed, transitions = _compute_transitions(curr, base)

    lines = []
    lines.append(f"Markers that changed category: {n_changed:,}")
    lines.append("")
    for (from_cat, to_cat), count in sorted(transitions.items(), key=lambda x: -x[1]):
        lines.append(f"  {from_cat:<32} → {to_cat:<32} {count:>8,}")
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    print(f"[benchmark] Loading manifest: {args.manifest}")
    main_df, chry_df, chr0_df = load_manifest(args.manifest)
    print(f"[benchmark]   Benchmarked={len(main_df):,}  Chr=Y={len(chry_df):,}  Chr=0={len(chr0_df):,}")

    assembly = detect_assembly(args.remapped, args.assembly)
    if args.assembly is None:
        print(f"[benchmark] Auto-detected assembly: {assembly}")

    print(f"[benchmark] Loading remapped: {args.remapped}")
    remapped_df = load_remapped(args.remapped, assembly)

    fasta = None
    if args.reference:
        import pysam
        print(f"[benchmark] Opening reference FASTA: {args.reference}")
        fasta = pysam.FastaFile(args.reference)

    try:
        print("[benchmark] Classifying markers...")
        result_df = compare_all(main_df, remapped_df, fasta=fasta, flank_len=args.flank_len)

        # Classify chrY and chr0 rows for their side TSVs (no explanatory layer — out of scope)
        chry_result = compare_all(chry_df, remapped_df)
        chry_result["result"] = "chrY"
        chr0_result = compare_all(chr0_df, remapped_df)
        chr0_result["result"] = "chr0"
    finally:
        if fasta is not None:
            fasta.close()

    # Merge QC filter trace if provided (benchmarked set only, per user scope).
    if args.traced:
        print(f"[benchmark] Loading QC trace: {args.traced}")
        trace_df = load_traced_why_filtered(args.traced, assembly)
        result_df = result_df.merge(trace_df, on="Name", how="left")
        result_df["why_filtered"] = result_df["why_filtered"].fillna("")

    # Write TSVs
    main_tsv = os.path.join(args.output_dir, f"benchmark_{ts}.tsv")
    chry_tsv = os.path.join(args.output_dir, f"benchmark_{ts}_chrY.tsv")
    chr0_tsv = os.path.join(args.output_dir, f"benchmark_{ts}_chr0.tsv")
    write_tsv(result_df,   main_tsv)
    write_tsv(chry_result, chry_tsv)
    write_tsv(chr0_result, chr0_tsv)
    print(f"[benchmark] TSVs written to: {args.output_dir}")

    # Load baseline if provided
    baseline_df = None
    if args.baseline:
        print(f"[benchmark] Loading baseline: {args.baseline}")
        baseline_df = load_baseline(args.baseline)

    # Write report
    report_path = os.path.join(args.output_dir, f"benchmark_{ts}_report.txt")
    write_report(
        result_df, chry_result, chr0_result,
        report_path,
        assembly=assembly,
        manifest_path=args.manifest,
        remapped_path=args.remapped,
        baseline_df=baseline_df,
    )
    print(f"[benchmark] Report: {report_path}")

    # Write standalone diff file if baseline was provided
    if baseline_df is not None:
        diff_path = os.path.join(args.output_dir, f"benchmark_{ts}_diff.txt")
        write_diff(result_df, baseline_df, diff_path)
        print(f"[benchmark] Diff:   {diff_path}")

    # Print report to stdout
    print()
    with open(report_path) as f:
        print(f.read())


if __name__ == "__main__":
    main()
