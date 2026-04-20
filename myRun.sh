srun --account=publicgrp -p low -t 01-0:00:00 -c 64 -n 1 -N 1 --mem=50g --pty bash

git clone git@github.com:drtamermansour/Equine80select_remapper.git
cd Equine80select_remapper
work_dir=$(pwd)


# Create the 'remap' conda environment with all dependencies
bash install.sh
conda activate remap

# Get the input files
parentageDir=$HOME/Horse_parentage_SNPs

# 1. Get the input manifest
mkdir -p manifests
cp $parentageDir/backup_original/Equine80select_24_20067593_B1.csv manifests/.
origManifest="$work_dir"/manifests/Equine80select_24_20067593_B1.csv
header_line=$(grep -n "^IlmnID" manifests/Equine80select_24_20067593_B1.csv | cut -d":" -f1)
end_line=$(grep -n "^\[Controls]" manifests/Equine80select_24_20067593_B1.csv | cut -d":" -f1)
nrows=$((end_line - header_line - 1)); echo $nrows ## 81974

## 2. get the reference genomes
## equCab3:
mkdir -p "$work_dir"/genomes/equCab3/download && cd "$work_dir"/genomes/equCab3/download
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.gz' -O equCab3.fa.gz
#gunzip equCab3.fa.gz
ln -s $parentageDir/equCab3/download/equCab3.fa .
cd "$work_dir"/genomes/equCab3
sed 's/>chr/>/' download/equCab3.fa > equCab3_genome.fa
equCab3_ref="$work_dir"/genomes/equCab3/equCab3_genome.fa
cd "$work_dir"


## find likely alternative haplotypes among unplaced scafolds
python scripts/scaffold_haplotype_analyzer.py \
	-r genomes/equCab3/equCab3_genome.fa \
	-o genomes/equCab3_scaffold_haplotype_analysis

python scripts/filter_scaffold_haplotypes.py \
	-i genomes/equCab3_scaffold_haplotype_analysis/scaffold_summary.tsv \
	-o genomes/equCab3_scaffold_haplotype_analysis/alt_haplotype_candidates.tsv \
	--min-identity 90 --min-query-cov 80 --max-span-ratio 3 --min-mapq 30 --max-blocks 5    

python scripts/exclude_alt_haplotypes.py \
    --scaffolds genomes/equCab3_scaffold_haplotype_analysis/alt_haplotype_candidates.tsv \
    --reference genomes/equCab3/equCab3_genome.fa \
    --output-dir genomes/equCab3_cleaned/   ## Removing 613 sequences → 4,088 retained ## cleaned FASTA: equCab3_cleaned/equCab3_genome_no_alt_haplotypes.fa


## This script is wrapper for "scripts/remap_manifest.py"
## Detailed pseudo-code for remap_manifest.py can be found in scripts/remap_manifest_psCode.txt
## scripts/remap_manifest.py add these columns to the manifest:
## Chr_EquCab3: chr on equCab3 based on 'TopGenomicSeq' alignment
## MapInfo_EquCab3: bp position on equCab3 (based primarily on Probe alignment; Fallback to TopGenomicSeq CIGAR.)
## Strand_EquCab3: SAM Flag from the 'TopGenomicSeq' alignment.
## Ref_EquCab3 & Alt_EquCab3: chosen from alleleA and alleleB obtained from 'TopGenomicSeq' e.g., "AGCT[A/G]TCGA"
## MAPQ_TopGenomicSeq: Mapping Quality score directly from the minimap2 alignment of the winning TopGenomicSeq candidate.
## MAPQ_Probe: The Mapping Quality score of the selected probe alignment. If no valid probe overlap was found (fallback used), this is set to 0.
#bash run_pipeline.sh \
#    --manifest manifests/Equine80select_24_20067593_B1.csv \
#    --reference genomes/equCab3/equCab3_genome.fa \
#    --assembly equCab3 \
#    --threads 64 \
#    --keep-temp \
#    --min-mapq-topseq 1 \
#    --resume \
#    --output-dir results_E80selv1_to_equCab3/

output_dir="results_E80selv2_to_equCab3noAlt_genDiv"

bash run_pipeline.sh \
    --manifest manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --reference genomes/equCab3_cleaned/equCab3_genome_no_alt_haplotypes.fa \
    --assembly equCab3noAlt \
    --threads 64 \
    --keep-temp \
    --resume \
    --output-dir "$output_dir" \
    --min-anchor topseq \
    --tie-policy avoid_scaffolds \
    --min-refalt-confidence low \
    --min-mapq-topseq 1
#    --preset permissive


## Running the Tests
pytest tests/ -v --manifest manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv --results-dir "$output_dir"


## Benchmarking Remapping Accuracy
python scripts/benchmark_compare.py \
    --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  "$output_dir"/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3noAlt.csv \
    --reference genomes/equCab3_cleaned/equCab3_genome_no_alt_haplotypes.fa \
    --traced "$output_dir"/qc/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3noAlt_traced.csv \
    --output-dir "$output_dir"/benchmark/

## compare accuracy of different coordinate systems 
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  "$output_dir"/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3noAlt.csv \
    --output-dir "$output_dir"/benchmark/



