srun --account=publicgrp -p low -t 01-0:0c0:00 -c 64 -n 1 -N 1 --mem=50g --pty bash

git clone git@github.com:drtamermansour/Equine80select_remapper.git
cd Equine80select_remapper
work_dir=$(pwd)
mkdir -p manifests
mkdir -p genomes

# Create the 'remap' conda environment with all dependencies
bash install.sh
conda activate remap

# Get the input files


## Benchmark aganist Illumina manifests as a gold standarad benchmark
## A) equCab3
# 1. Get the input manifest
module load rclone ## Loading rclone/1.65.1
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/Miscellaneous documents_standardbred/Equine80select_v2_1_HTS_20143333_B1_UCD.csv" --drive-shared-with-me manifests/.
origManifest="$work_dir"/manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv
header_line=$(grep -n "^IlmnID" "$origManifest" | cut -d":" -f1)
end_line=$(grep -n "^\[Controls]" "$origManifest" | cut -d":" -f1)
rows=$((end_line - header_line - 1)); echo $nrows ## 84319

## 2. get the reference genomes
## equCab3:
mkdir -p "$work_dir"/genomes/equCab3/download && cd "$work_dir"/genomes/equCab3/download
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.gz' -O equCab3.fa.gz
#gunzip equCab3.fa.gz
parentageDir=$HOME/Horse_parentage_SNPs
ln -s $parentageDir/genomes/equCab3/download/equCab3.fa .
cd "$work_dir"/genomes/equCab3
sed 's/>chr/>/' download/equCab3.fa > equCab3_genome.fa
equCab3_ref="$work_dir"/genomes/equCab3/equCab3_genome.fa
cd "$work_dir"



## Analysis scenario 1: Baseline analysis 
# Run the pipeline using the raw reference (no clean up of alt haplotypes), and apply no filters (Except those enforced) 
bash run_pipeline.sh \
    --manifest manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --reference genomes/equCab3/equCab3_genome.fa \
    --assembly equCab3 \
    --threads 64 \
    --keep-temp \
    --resume \
    --output-dir results_E80selv2_to_equCab3_noFilters/ \
    --preset permissive


## Running the Tests
pytest tests/ -v --manifest manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv --results-dir results_E80selv2_to_equCab3_noFilters


## Benchmarking Remapping Accuracy
python scripts/benchmark_compare.py \
    --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3_noFilters/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
    --reference genomes/equCab3/equCab3_genome.fa \
    --traced results_E80selv2_to_equCab3_noFilters/qc/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3_traced.csv \
    --output-dir results_E80selv2_to_equCab3_noFilters/benchmark/


## compare accuracy of different coordinate systems 
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3_noFilters/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3.csv \
    --output-dir results_E80selv2_to_equCab3_noFilters/benchmark/



#####################################################################################################
## Analysis scenario 2: Using cleaned reference 
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

bash run_pipeline.sh \
    --manifest manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --reference genomes/equCab3_cleaned/equCab3_genome_no_alt_haplotypes.fa \
    --assembly equCab3noAlt \
    --threads 64 \
    --keep-temp \
    --resume \
    --output-dir results_E80selv2_to_equCab3noAlt_noFilters/ \
    --preset permissive


## Running the Tests
pytest tests/ -v --manifest manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv --results-dir results_E80selv2_to_equCab3noAlt_noFilters


## Benchmarking Remapping Accuracy
python scripts/benchmark_compare.py \
    --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3noAlt_noFilters/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3noAlt.csv \
    --reference genomes/equCab3_cleaned/equCab3_genome_no_alt_haplotypes.fa \
    --traced results_E80selv2_to_equCab3noAlt_noFilters/qc/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3noAlt_traced.csv \
    --output-dir results_E80selv2_to_equCab3noAlt_noFilters/benchmark/

## compare accuracy of different coordinate systems 
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --remapped  results_E80selv2_to_equCab3noAlt_noFilters/remapping/Equine80select_v2_1_HTS_20143333_B1_UCD_remapped_equCab3noAlt.csv \
    --output-dir results_E80selv2_to_equCab3noAlt_noFilters/benchmark/

########################################
## B) canFam3
# 1. Get the input manifest
## Infinium CanineHTS BeadChip
## https://support.illumina.com/downloads/caninehts_product_files.html


## 2. get the reference genomes
## equCab3:
mkdir -p "$work_dir"/genomes/canFam3/download && cd "$work_dir"/genomes/canFam3/download
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz' -O canFam3.fa.gz
#gunzip canFam3.fa.gz
canineDir=$HOME/MAF
ln -s $canineDir/refGenome/canFam3.fa .
cd "$work_dir"/genomes/canFam3
sed 's/>chr/>/' download/canFam3.fa > canFam3_genome.fa
canFam3_ref="$work_dir"/genomes/canFam3/canFam3_genome.fa
cd "$work_dir"


## Analysis scenario 1: Baseline analysis 
# Run the pipeline using the raw reference (no clean up of alt haplotypes), and apply no filters (Except those enforced) 
man_prefix="CanineHTS-24_20095584_B1"
ref="canFam3"

bash run_pipeline.sh \
    --manifest manifests/CanineHTS-24_20095584_B1.csv \
    --reference genomes/canFam3/canFam3_genome.fa \
    --assembly canFam3 \
    --threads 64 \
    --keep-temp \
    --resume \
    --output-dir results_CanineHTS24_to_canFam3_noFilters/ \
    --preset permissive

## Running the Tests
pytest tests/ -v --manifest manifests/CanineHTS-24_20095584_B1.csv --results-dir results_CanineHTS24_to_canFam3_noFilters


## Benchmarking Remapping Accuracy
python scripts/benchmark_compare.py \
    --manifest  manifests/CanineHTS-24_20095584_B1.csv \
    --remapped  results_CanineHTS24_to_canFam3_noFilters/remapping/CanineHTS-24_20095584_B1_remapped_canFam3.csv \
    --reference genomes/canFam3/canFam3_genome.fa \
    --traced results_CanineHTS24_to_canFam3_noFilters/qc/CanineHTS-24_20095584_B1_remapped_canFam3_traced.csv \
    --output-dir results_CanineHTS24_to_canFam3_noFilters/benchmark/


## compare accuracy of different coordinate systems 
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  manifests/CanineHTS-24_20095584_B1.csv \
    --remapped  results_CanineHTS24_to_canFam3_noFilters/remapping/CanineHTS-24_20095584_B1_remapped_canFam3.csv \
    --output-dir results_CanineHTS24_to_canFam3_noFilters/benchmark/

#############################################
## susScr11
# 1. Get the input manifest
## PorcineSNP60 v3.0 BeadChip
## https://support.illumina.com/downloads/porcinesnp60-v3_product_files.html

## 2. get the reference genomes
## susScr11:
mkdir -p "$work_dir"/genomes/susScr11/download && cd "$work_dir"/genomes/susScr11/download
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/susScr11/bigZips/susScr11.fa.gz' -O susScr11.fa.gz
gunzip susScr11.fa.gz
cd "$work_dir"/genomes/susScr11
sed 's/>chr/>/' download/susScr11.fa > susScr11_genome.fa
equCab3_ref="$work_dir"/genomes/susScr11/susScr11_genome.fa
cd "$work_dir"

## Analysis scenario 1: Baseline analysis 
# Run the pipeline using the raw reference (no clean up of alt haplotypes), and apply no filters (Except those enforced) 
man_prefix="PorcineSNP80v1_HTS_20033000_A2"
ref="susScr11"
out="PorcineSNP60_to_susScr11"

bash run_pipeline.sh \
    --manifest manifests/$man_prefix.csv \
    --reference genomes/${ref}/${ref}_genome.fa \
    --assembly ${ref} \
    --threads 64 \
    --keep-temp \
    --resume \
    --output-dir results_${out}_noFilters/ \
    --preset permissive

## Running the Tests
pytest tests/ -v --manifest manifests/$man_prefix.csv --results-dir results_${out}_noFilters


## Benchmarking Remapping Accuracy
python scripts/benchmark_compare.py \
    --manifest  manifests/$man_prefix.csv \
    --remapped  results_${out}_noFilters/remapping/${man_prefix}_remapped_${ref}.csv \
    --reference genomes/${ref}/${ref}_genome.fa \
    --traced results_${out}_noFilters/qc/${man_prefix}_remapped_${ref}_traced.csv \
    --output-dir results_${out}_noFilters/benchmark/


## compare accuracy of different coordinate systems 
python scripts/benchmark_cigar_vs_probe.py \
    --manifest  manifests/$man_prefix.csv \
    --remapped  results_${out}_noFilters/remapping/${man_prefix}_remapped_${ref}.csv \
    --output-dir results_${out}_noFilters/benchmark/

##################################################
## Benchmark aganist liftover/CrossMap
## Already created genomes/equCab3_cleaned/equCab3_genome_no_alt_haplotypes.fa

bash scripts/benchmark/run_benchmark_vs_liftover.sh \
    --v1         manifests/Equine80select_24_20067593_B1.csv \
    --v2         manifests/Equine80select_v2_1_HTS_20143333_B1_UCD.csv \
    --reference  genomes/equCab3/equCab3_genome.fa \
    --output-dir results_liftover_benchmark/ \
    --threads    64 \
    --resume                                              



