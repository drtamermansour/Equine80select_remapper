"""
exclude_alt_haplotypes.py — Build a cleaned reference FASTA by removing scaffolds
identified as alternative haplotypes by filter_scaffold_haplotypes.py.

Reads the scaffold list from the filter output (TSV with scaffold_id as first column),
removes those sequences from the reference FASTA, writes a new indexed FASTA to the
output directory, and reports counts.

Usage:
  python scripts/exclude_alt_haplotypes.py \\
      --scaffolds remap_assessment/scaffold_haplotype_analysis/alt_haplotype_candidates.tsv \\
      --reference equCab3/equCab3_genome.fa \\
      --output-dir equCab3_cleaned/

Outputs (in --output-dir):
  {stem}_no_alt_haplotypes.fa       new FASTA (scaffolds excluded)
  {stem}_no_alt_haplotypes.fa.fai   samtools index
  exclusion_report.txt              summary of excluded / retained sequences
"""

import argparse
import os
import subprocess
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Remove alt-haplotype scaffolds from a reference FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--scaffolds",   required=True,
                   help="TSV from filter_scaffold_haplotypes.py; scaffold_id must be the first column")
    p.add_argument("--reference",   required=True,
                   help="Reference FASTA to filter (must have an .fai beside it)")
    p.add_argument("--output-dir",  required=True,
                   help="Directory for the cleaned FASTA and index (created if absent)")
    return p.parse_args()


def read_scaffold_ids(tsv_path):
    """
    Read the scaffold_id column from the filter_scaffold_haplotypes.py output TSV.
    Skips the header row. The first column is always scaffold_id.
    Returns a set of scaffold IDs to exclude.
    """
    ids = set()
    with open(tsv_path) as fh:
        header = fh.readline().strip().split("\t")
        if not header or header[0] != "scaffold_id":
            raise ValueError(
                f"Expected first column to be 'scaffold_id', got {header[0]!r}.\n"
                f"Ensure --scaffolds points to filter_scaffold_haplotypes.py output."
            )
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ids.add(line.split("\t")[0])
    return ids


def read_fai(fai_path):
    """Return list of all sequence names in the FAI, in FAI order."""
    names = []
    with open(fai_path) as fh:
        for line in fh:
            cols = line.split("\t")
            if cols:
                names.append(cols[0])
    return names


def extract_sequences(seq_names, ref_fa, out_fa):
    """
    Extract named sequences from ref_fa into out_fa using samtools faidx.
    Sequences are written in the order given by seq_names.
    """
    cmd = ["samtools", "faidx", ref_fa] + seq_names
    with open(out_fa, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit(f"[exclude_alt] ERROR: samtools faidx failed:\n{result.stderr}")


def index_fasta(fa_path):
    """Run samtools faidx to build the .fai index for the new FASTA."""
    result = subprocess.run(
        ["samtools", "faidx", fa_path],
        stderr=subprocess.PIPE, text=True,
    )
    if result.returncode != 0:
        sys.exit(f"[exclude_alt] ERROR: samtools faidx indexing failed:\n{result.stderr}")


def write_report(report_path, ref_fa, scaffolds_tsv, out_fa,
                 all_names, excluded, retained):
    lines = [
        "Exclusion Report — exclude_alt_haplotypes.py",
        f"Reference:        {ref_fa}",
        f"Scaffold list:    {scaffolds_tsv}",
        f"Output FASTA:     {out_fa}",
        "-" * 60,
        f"Total sequences in reference:  {len(all_names):>8,}",
        f"Scaffolds in exclusion list:   {len(excluded):>8,}",
        f"  of which found in reference: {len(excluded & set(all_names)):>8,}",
        f"  not found (already absent):  {len(excluded - set(all_names)):>8,}",
        f"Sequences retained in output:  {len(retained):>8,}",
        "",
    ]
    if excluded - set(all_names):
        lines.append("Scaffold IDs requested for exclusion but not in reference:")
        for name in sorted(excluded - set(all_names)):
            lines.append(f"  {name}")
        lines.append("")
    with open(report_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    args = parse_args()

    # Validate inputs
    fai_path = args.reference + ".fai"
    if not os.path.exists(args.reference):
        sys.exit(f"[exclude_alt] ERROR: Reference not found: {args.reference}")
    if not os.path.exists(fai_path):
        sys.exit(
            f"[exclude_alt] ERROR: FAI index not found: {fai_path}\n"
            f"Run: samtools faidx {args.reference}"
        )
    if not os.path.exists(args.scaffolds):
        sys.exit(f"[exclude_alt] ERROR: Scaffold list not found: {args.scaffolds}")

    os.makedirs(args.output_dir, exist_ok=True)

    # Derive output filename from reference stem
    ref_stem = os.path.splitext(os.path.basename(args.reference))[0]
    # strip .fa / .fasta extension if double-extension (e.g. genome.fa.gz handled upstream)
    out_fa = os.path.join(args.output_dir, f"{ref_stem}_no_alt_haplotypes.fa")
    report_path = os.path.join(args.output_dir, "exclusion_report.txt")

    # Read scaffold IDs to exclude
    print(f"[exclude_alt] Reading scaffold list: {args.scaffolds}")
    excluded = read_scaffold_ids(args.scaffolds)
    print(f"[exclude_alt]   {len(excluded):,} scaffolds to exclude")

    # Read all sequence names from FAI
    all_names = read_fai(fai_path)
    print(f"[exclude_alt] Reference has {len(all_names):,} sequences")

    # Compute retained set (preserve FAI order)
    retained = [name for name in all_names if name not in excluded]
    n_removed = len(all_names) - len(retained)
    not_found = excluded - set(all_names)

    if not_found:
        print(f"[exclude_alt] WARNING: {len(not_found):,} scaffold IDs not found in reference "
              f"(may have been already absent).")
    print(f"[exclude_alt] Removing {n_removed:,} sequences → {len(retained):,} retained")

    if not retained:
        sys.exit("[exclude_alt] ERROR: All sequences would be removed. Check --scaffolds input.")

    # Extract retained sequences
    print(f"[exclude_alt] Writing cleaned FASTA: {out_fa}")
    extract_sequences(retained, args.reference, out_fa)

    # Index the new FASTA
    print(f"[exclude_alt] Indexing: {out_fa}.fai")
    index_fasta(out_fa)

    # Write report
    write_report(report_path, args.reference, args.scaffolds, out_fa,
                 all_names, excluded, retained)
    print(f"[exclude_alt] Report: {report_path}")
    print(f"[exclude_alt] Done. Cleaned reference: {out_fa}")


if __name__ == "__main__":
    main()
