# CLAUDE.md — Developer Quick Reference

Compact orientation for future Claude Code sessions. User-facing docs live in
`README.md` and `docs/`; this file is for fast retrieval of invariants,
pitfalls, file paths, and test commands.

## What this is

Pipeline that remaps Illumina genotyping arrays between reference assemblies.
Built for **Equine80select v1 (EquCab2) → EquCab3** (~82 k SNPs).

Headline output: `qc/matchingSNPs_binary_consistantMapping.{assembly}_map`
(8-column TSV, no header — see `docs/output_formats.md`).

## File map

```
run_pipeline.sh                    orchestrator: faidx → remap_manifest.py → qc_filter.py
submit_slurm.sh                    SLURM wrapper
install.sh                         conda env 'remap' (minimap2, samtools, bcftools, pysam, pandas)
myRun.sh                           historical record of the v1→equCab3 invocation

scripts/_strand_utils.py           strand_normalize(allele, strand) + complement(seq)
scripts/remap_manifest.py          dual alignment + 3-dim output
scripts/qc_filter.py               11-stage QC cascade + WhyFiltered trace
scripts/benchmark_compare.py       ground-truth benchmark + verdict + QC impact
scripts/benchmark_cigar_vs_probe.py  3-way coord comparison (probe_cigar / topseq_cigar / final)
scripts/scaffold_haplotype_*.py    pre-pipeline alt-haplotype removal
scripts/exclude_alt_haplotypes.py  builds cleaned FASTA

tests/                             pytest; conftest takes --results-dir for integration
docs/                              user-facing guides; *_simple.md is the source mermaid
```

## 3-dimension output framework (single source of truth)

Every marker gets exactly three orthogonal columns (after `-a equCab3`):

- `anchor_equCab3` ∈ {`topseq_n_probe`, `topseq_only`, `probe_only`, `N/A`}
- `tie_equCab3` ∈ {`unique`, `AS_resolved`, `dAS_resolved`, `NM_resolved`,
  `CoordDelta_resolved`, `scaffold_resolved`, `locus_unresolved`, `N/A`}
- `RefAltMethodAgreement_equCab3` ∈ {`NM_match`, `NM_validated`, `NM_N/A`,
  `NM_tied`, `NM_only`, `NM_unmatch`, `NM_corrected`, `NM_mismatch`,
  `refalt_unresolved`, `N/A`}

**Chr=0 rule:** marker is excluded from downstream iff `anchor=N/A` OR
`tie=locus_unresolved` OR `RefAltMethodAgreement=refalt_unresolved` (or
`NM_mismatch` for indels via Stage 2 of QC).

## QC cascade (11 stages, all configurable)

In `qc_filter.py:WHY_FILTERED_LABELS`. Order matters because `WhyFiltered`
records the **first** rejection per marker:

1. `stage_1_failed_markers` — `Strand=N/A` (always on)
2. `stage_2_design_conflict` — Ref ≠ genome ref / `NM_mismatch` (always on)
3. `stage_3_coordinate_role` — `--coordinate-role` (`High`/`Moderate`/`Low`)
4. `stage_4_tie_label` — `--tie-label` (`unique`/`resolved`/`avoid_scaffolds`)
5. `stage_5_refalt_conf` — `--refalt-conf` (`High`/`Moderate`/`Low`)
6. `stage_6_mapq_topseq` — `--mapq-topseq 30` (probe_only NaN-exempt)
7. `stage_7_mapq_probe` — `--mapq-probe 0` disabled (topseq_only NaN-exempt)
8. `stage_8_coord_delta` — `--coord-delta -1` disabled
9. `stage_9_indel_excluded` — off by default; `--keep-indels` to include
10. `stage_10_polymorphic` — off by default; `--keep-polymorphic` to include
11. `stage_11_ambiguous_snp` — off by default; `--keep-ambiguous` to include

`tie-label resolved` accepts all `*_resolved` **except** `scaffold_resolved`
(see `_TIE_RESOLVED` in `qc_filter.py`).

## Coordinate selection (in `run_remapping`)

Two CIGAR-derived coordinates per `topseq_n_probe` marker:

- `CoordProbe_{a}` — from probe alignment's CIGAR walk.
- `Coord_TopSeqCIGAR_{a}` — from TopGenomicSeq alignment's CIGAR walk.
- `CoordDelta_{a}` = `|CoordProbe − Coord_TopSeqCIGAR|`; `−1` if CIGAR unavailable.
- `CoordSource_{a}` ∈ {`"probe_cigar"`, `"topseq_cigar"`, `"N/A"`} — which
  one ended up in `MapInfo_{a}`.

Rule: indels OR `CoordDelta ≥ 2` → use `topseq_cigar`. Otherwise → `probe_cigar`.
`CoordDelta = -1` (soft clip) → `probe_cigar`. Empirical accuracy:
probe_cigar 98.0%, topseq_cigar 98.6%, final 98.7%.

## Critical invariants & pitfalls

1. **Pipeline-internal code must NOT read `RefStrand`/`SourceStrand`/`IlmnStrand` as ground truth.** They are Illumina-design conventions that don't reliably predict reference-strand orientation. Only `benchmark_compare.py` uses `RefStrand` (as benchmark ground truth). Strand handling inside `remap_manifest.py` is sequence-derived (`probe_topseq_orientation` + 21-mer fallback).

2. **`MAPQ_Probe = NaN` ≠ MAPQ was zero.** NaN means `topseq_only` (no probe alignment used). MAPQ filters must do `.notna() & (col < threshold)` — see `apply_probe_mapq_filter`. Same for `MAPQ_TopGenomicSeq = NaN` for `probe_only` markers.

3. **`CoordDelta = -1` is a sentinel**, not a negative distance. `topseq_only` and `probe_only` carry it because there's no probe-vs-CIGAR pair to compare. Stage 8 `--coord-delta N` correctly leaves them through (`-1 > N` is false for N ≥ 0).

4. **Probe alignment strand ≠ TopSeq alignment strand.** They are independent SAM-FLAG-derived values. Mixing them caused the original ±51 bp bug. `get_probe_coordinate` uses **probe** strand; `parse_cigar_to_ref_pos` uses **TopSeq** strand for `target_idx`.

5. **`Ref_{a}` / `Alt_{a}` are stored in the TopSeq alignment-strand orientation.** `qc_filter.py` strand-normalises them to `+` strand for VCF/BIM output via `strand_normalize` from `_strand_utils`.

6. **Deletion minus-strand correction** in `run_remapping` (after `determine_ref_alt_v2`): when `len(ref) > len(alt)` and TopSeq strand is `-`, subtract `len(ref) - len(alt)` from `final_pos` and `c_pos` if `coord_source == "probe_cigar"`. Don't touch this — it's correct.

7. **Resume after step 3 failure** — `bash run_pipeline.sh ... --resume --temp-dir <existing>` skips the multi-hour minimap2 step. Use `--keep-temp` to preserve SAM files first.

8. **`Chr_{a}` is always a string.** Load with `dtype={col_chr: str}`. Unplaced contigs like `Un_NW_*` are valid values.

9. **`bcftools` must be on PATH** for `qc_filter.py extract_ref_alleles`. `conda activate remap` provides it.

10. **`-a / --assembly` is required** in `run_pipeline.sh` (no FASTA-derived default — column names embed the assembly label and a wrong derivation propagates silently).

## Tests

```bash
conda activate remap
pytest tests/ -v                                         # 271 unit tests (~9 s)
pytest tests/test_benchmark_compare.py --results-dir <dir> -v   # +3 integration
```

When adding production code, mirror the existing TDD pattern (RED → GREEN
in `tests/test_<module>.py`). Tests use `_strand_utils.strand_normalize`
directly via `scripts/` import (see `tests/conftest.py` for the `sys.path`
manipulation).

## Common dev tasks

- **Adding a QC filter stage** — extend `WHY_FILTERED_LABELS` (in order),
  add the helper function (mirror `apply_exclude_indels_filter`), wire into
  `run_qc` after the previous `_tag_removed` call, add CLI flag in
  `parse_args`, plumb through `run_pipeline.sh`, add row to README's
  cascade table and the `docs/qc_filters.md` table.

- **Adding a benchmark verdict** — add the rule in
  `benchmark_compare.classify_explanatory`, update the verdict table in
  `docs/benchmarking.md`, add a unit test in `tests/test_benchmark_compare.py`.

- **Touching strand normalisation** — update only `scripts/_strand_utils.py`.
  All consumers (`remap_manifest.py`, `qc_filter.py`, `benchmark_compare.py`)
  import from it. Don't reintroduce a hand-rolled `_RC_TABLE` or
  `COMPLEMENT` dict in any consumer.

## Things to ignore

- `temp/` and `ignore/` — historical scripts and superseded plans.
- Any file under `ignore/superpowers/` — past brainstorm/plan artifacts.
- `my-session.txt` — historical session log.
- `backup_original/` — original Illumina manifests; do not modify.
