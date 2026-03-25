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
mkdir -p backup_original
cp $parentageDir/backup_original/Equine80select_24_20067593_B1.csv backup_original/.
origManifest="$work_dir"/backup_original/Equine80select_24_20067593_B1.csv
header_line=$(grep -n "^IlmnID" backup_original/Equine80select_24_20067593_B1.csv | cut -d":" -f1)
end_line=$(grep -n "^\[Controls]" backup_original/Equine80select_24_20067593_B1.csv | cut -d":" -f1)
nrows=$((end_line - header_line - 1)); echo $nrows ## 81974

## 2. get the reference genomes
## equCab3:
mkdir -p "$work_dir"/equCab3/download && cd "$work_dir"/equCab3/download
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.gz' -O equCab3.fa.gz
#gunzip equCab3.fa.gz
ln -s $parentageDir/equCab3/download/equCab3.fa .
cd "$work_dir"/equCab3
sed 's/>chr/>/' download/equCab3.fa > equCab3_genome.fa
equCab3_ref="$work_dir"/equCab3/equCab3_genome.fa
cd "$work_dir"



## This script is wrapper for "scripts/remap_manifest.py"
## Detailed pseudo-code for remap_manifest.py can be found in scripts/remap_manifest_psCode.txt
## scripts/remap_manifest.py add these columns to the manifest:
## Chr_EquCab3: chr on equCab3 based on 'TopGenomicSeq' alignment
## MapInfo_EquCab3: bp position on equCab3 (based primarily on Probe alignment; Fallback to TopGenomicSeq CIGAR.)
## Strand_EquCab3: SAM Flag from the 'TopGenomicSeq' alignment.
## Ref_EquCab3 & Alt_EquCab3: chosen from alleleA and alleleB obtained from 'TopGenomicSeq' e.g., "AGCT[A/G]TCGA"
## MAPQ_TopGenomicSeq: Mapping Quality score directly from the minimap2 alignment of the winning TopGenomicSeq candidate.
## MAPQ_Probe: The Mapping Quality score of the selected probe alignment. If no valid probe overlap was found (fallback used), this is set to 0.
bash run_pipeline.sh \
    --manifest backup_original/Equine80select_24_20067593_B1.csv \
    --reference equCab3/equCab3_genome.fa \
    --assembly equCab3 \
    --threads 64 \
    --keep-temp \
    --mapq-topseq 1 \
    --resume \
    --output-dir results/



