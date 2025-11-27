##!/bin/bash
#SBATCH --job-name=remap_equCab3
#SBATCH --output=remap_equCab3_%j.out
#SBATCH --error=remap_equCab3_%j.err
#SBATCH --partition=bmm
#SBATCH --time=01-0:00:00
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --ntasks=1 
srun -p bmm -t 01-0:00:00 -c 64 -n 1 -N 1 --mem=50g --pty bash

## create a new conda environment for genomic related GWAS tools
mamba create -n remap -c bioconda minimap2 
conda activate remap
mamba install -c bioconda pysam samtools bcftools
mamba install -c conda-forge pandas

mkdir -p equCab3_remap && cd equCab3_remap
work_dir=$(pwd)
parentageDir=$work_dir/../Horse_parentage_SNPs

mkdir -p backup_original
cp $parentageDir/backup_original/Equine80select_24_20067593_B1.csv backup_original/.
origManifest="$work_dir"/backup_original/Equine80select_24_20067593_B1.csv
header_line=$(grep -n "^IlmnID" backup_original/Equine80select_24_20067593_B1.csv | cut -d":" -f1)
end_line=$(grep -n "^\[Controls]" backup_original/Equine80select_24_20067593_B1.csv | cut -d":" -f1)
nrows=$((end_line - header_line - 1)); echo $nrows ## 81974

## get the reference genomes
## equCab2:
mkdir -p "$work_dir"/equCab2/download && cd "$work_dir"/equCab2/download
#wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab2/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
#tar xvzf chromFa.tar.gz
cd "$work_dir"/equCab2
#cat $(ls -1v download/*.fa) > equCab2_genome.fa
#sed -i 's/>chr/>/' equCab2_genome.fa
ln -s $parentageDir/equCab2/equCab2_genome.fa .
equCab2_ref="$work_dir"/equCab2/equCab2_genome.fa
cat $equCab2_ref | grep -v "^#" | awk -F">" '{if(!$1){if(NR!=1)print "##contig=<ID="ref",length="reflen">";ref=$2;reflen=0}else reflen+=length($1)}END{print "##contig=<ID="ref",length="reflen">";}' > vcf_contigs.txt
equCab2_vcfContigs="$work_dir"/equCab2/vcf_contigs.txt

## equCab3:
mkdir -p "$work_dir"/equCab3/download && cd "$work_dir"/equCab3/download
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.gz' -O equCab3.fa.gz
#gunzip equCab3.fa.gz
ln -s $parentageDir/equCab3/download/equCab3.fa .
cd "$work_dir"/equCab3
sed 's/>chr/>/' download/equCab3.fa > equCab3_genome.fa
equCab3_ref="$work_dir"/equCab3/equCab3_genome.fa
cat $equCab3_ref | grep -v "^#" | awk -F">" '{if(!$1){if(NR!=1)print "##contig=<ID="ref",length="reflen">";ref=$2;reflen=0}else reflen+=length($1)}END{print "##contig=<ID="ref",length="reflen">";}' > vcf_contigs.txt
equCab3_vcfContigs="$work_dir"/equCab3/vcf_contigs.txt

# create samtools faidx index
samtools faidx $equCab3_ref
##################################

# Run the remapping steps
cd "$work_dir"

#mkdir tmpX && mv *.txt *remapped* temp_* equCab3.mismatches tmpX/.
python scripts/remap_manifest.py \
    -i $origManifest \
    -r $equCab3_ref \
    -o Equine80select_24_20067593_B1_remapped.csv

#####################
## Explore the manifest header and format
mkdir -p temp_explore
grep IlmnID $origManifest | tr ',' '\n' > temp_explore/1.txt
grep IlmnID Equine80select_24_20067593_B1_remapped.csv | tr ',' '\n' | awk '{print NR}' > temp_explore/2.txt
grep IlmnID Equine80select_24_20067593_B1_remapped.csv | tr ',' '\n' > temp_explore/3.txt
grep -A1 IlmnID Equine80select_24_20067593_B1_remapped.csv | tail -n1 | tr ',' '\n' > temp_explore/4.txt
paste temp_explore/1.txt temp_explore/2.txt temp_explore/3.txt temp_explore/4.txt


#####################
## Benchmark aganist the already known EquCab3 markers 
mkdir -p remap_assessment
## Summary statistics of input EquCab2 and EquCab3 markers of the manifest
cat Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $9}' | sort | uniq -c
#  75952 2.0
#   6022 3.0

cat Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $20}' | sort | uniq -c | awk '{if($2==0) print}'  ## 29 with 0 value (i.e. failed to remap)
cat Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $20}' | grep Un | wc -l ## 644 assigned to chrUn

#cat Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{if($9==3)print $4,$3"/"$16,$22}' | sort | uniq -c > remap_assessment/equCab3.metadata_counts
cat Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{if($9==3)print $2,$4,$10":"$11,$3"/"$16,$20":"$21,$22}' | awk 'BEGIN{FS=","}{if($3==$5)print}' | wc -l ## 5949 remapped correctly
cat Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{if($9==3)print $2,$4,$10":"$11,$3"/"$16,$20":"$21,$22}' | awk 'BEGIN{FS=","}{if($3!=$5)print}' > remap_assessment/equCab3.mismatches
cat remap_assessment/equCab3.mismatches | wc -l ## 73
cat remap_assessment/equCab3.mismatches | awk 'BEGIN{FS=OFS=","}{print $2,$4,$6}' | sort | uniq -c
#      6 [I/D],MINUS/PLUS,+
#     16 [I/D],PLUS/PLUS,+
#      1 [A/C],TOP/TOP,+
#      5 [A/G],TOP/BOT,+
#      5 [A/G],TOP/BOT,-
#      3 [A/G],TOP/TOP,+
#      3 [A/G],TOP/TOP,-
#      5 [A/T],TOP/BOT,-
#      1 [C/G],TOP/TOP,-
#      3 [G/C],BOT/TOP,+
#      2 [G/C],BOT/TOP,-
#      5 [T/C],BOT/BOT,+
#     10 [T/C],BOT/BOT,-
#      3 [T/C],BOT/TOP,-
#      1 [T/G],BOT/BOT,+
#      1 [T/G],BOT/BOT,-
#      1 [T/G],BOT/TOP,+
#      2 [T/G],BOT/TOP,-

#######################
## Extract discrepancies to separate files
head -n1 Equine80select_24_20067593_B1_remapped.csv > remap_assessment/remapped_discrepancies.csv
cat remap_assessment/equCab3.mismatches | cut -f1 -d, | grep -Fwf - Equine80select_24_20067593_B1_remapped.csv >> remap_assessment/remapped_discrepancies.csv
grep ^@ temp_probe.sam > remap_assessment/remapped_discrepancies_probe.sam
cat remap_assessment/equCab3.mismatches | cut -f1 -d, | grep -Fwf - temp_probe.sam >> remap_assessment/remapped_discrepancies_probe.sam
grep ^@ temp_topseq.sam > remap_assessment/remapped_discrepancies_topseq.sam
cat remap_assessment/equCab3.mismatches | cut -f1 -d, | grep -Fwf - temp_topseq.sam >> remap_assessment/remapped_discrepancies_topseq.sam

module load rclone
rclone -v --copy-links copy remapped_discrepancies.* remote_UCDavis_GoogleDr:temp/
##########################
## Generate histogram of MAPQ scores for all TopGenomicSeq alignments
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($25/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 Equine80select_24_20067593_B1_remapped.csv) > remap_assessment/MAPQ_TopGenomicSeq.histo

## Generate histogram of MAPQ scores for all probes alignments
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($26/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 Equine80select_24_20067593_B1_remapped.csv) > remap_assessment/MAPQ_Probe.histo

cat <(grep -v "^@" temp_topseq.sam | sed 's/_A\t/\t/;s/_B\t/\t/;') <(grep -v "^@" temp_probe.sam) | cut -f1,3 | sort | uniq -c | awk '{if($1!=3)print}' > inconsistant_alignments 
cat inconsistant_alignments | awk '{print $2}' | sort | uniq > inconsistent_probes
cat inconsistent_probes | grep -Fwf - Equine80select_24_20067593_B1_remapped.csv > inconsistent_remapped.csv

## histogram of MAPQ scores for TopGenomicSeq alignments from inconsistent_probes
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($25/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 inconsistent_remapped.csv) > remap_assessment/MAPQ_TopGenomicSeq_inconsistent_remapped.histo

## Generate histogram of MAPQ scores for probes with inconsistent alignments
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($26/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 inconsistent_remapped.csv) > remap_assessment/MAPQ_Probe_inconsistent_remapped.histo

##########################
## Assessemnt of the remapping results
#echo "SNP_ID,AlleleA,AlleleB,IlmnStrand,SourceStrand,Ref_Allele,Alt_Allele,strand" | tr ',' '\t' > remap_assessment/allele_assessment.tsv
#tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=",";OFS="\t"}{gsub(/\[|\]/, "", $4);gsub(/\//, OFS, $4);print $2,$4,$3,$16,$23,$24,$22}' >> remap_assessment/allele_assessment.tsv
#python scripts/validate_alleleMatch_v2.py --input remap_assessment/allele_assessment.tsv --output remap_assessment/allele_assessment_validated.tsv

# tranform the array manifest into VCF & generate ALT/REF map 
# Modified from 'manifest_to_ref.sh'
#cat Equine80select_24_20067593_B1_remapped.csv |\
#    awk 'BEGIN{FS=OFS=","} NR==1{print; next} {if ($9==3 && ($10 FS $11)==($20 FS $21)) print}' > remapped_markers.csv ## 5949
cat Equine80select_24_20067593_B1_remapped.csv |\
    awk 'BEGIN{FS=OFS=","} NR==1{print; next} {if ($20!=0) print}' > remapped_markers.csv ## 81945
manifest="remapped_markers.csv"
vcfContigs="$equCab3_vcfContigs"
ref="$equCab3_ref"
out="Equine80select_24_20067593_B1_remapped"
# create VCF template using the positions of the array SNPs
echo "##fileformat=VCFv4.3" > _pos.vcf
cat $vcfContigs >> _pos.vcf
echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr ' ' '\t' >> _pos.vcf
tail -n+2 $manifest | awk 'BEGIN{FS=",";OFS="\t"}{print $20,$21,$2,"N",".",".",".","."}' >> _pos.vcf

# Use bcftools to obtain the reference alleles of the array SNPs
bcftools norm -c ws -f $ref _pos.vcf 1> _ref.vcf 2> _ref.vcf.log    # total/modified/added:   81945/0/81945

# align the vcf carring reference alleles with the vcf carring positions (which already matches the manifest)
cp _ref.vcf _ref.vcf_temp
grep "^#" _ref.vcf_temp > _ref.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$4;next;}{$4=a[$3];print}' <(grep -v "^#" _ref.vcf_temp) <(grep -v "^#" _pos.vcf) >> _ref.vcf
rm _ref.vcf_temp

# Transform the array SNP alleles into reference and alternative alleles
grep -v "^#" _ref.vcf | awk 'BEGIN{FS=OFS="\t"}{print $4}' | tr 'tcga' 'TCGA' > _tmp.ref_alleles # single column of ref allele
tail -n+2 "$manifest" | awk 'BEGIN{FS=",";OFS="\t"}{print $20,$21,$2,$22,$23,$24}' > _tmp.mapped_alleles # the new alleles and their strand # 11  21962991  21962991_Curly_f_ilmndup1   +   G   A
tail -n+2 "$manifest" | awk 'BEGIN{FS=",";OFS="\t"}{print $23,$24}' | tr 'TCGA' 'AGCT' > _tmp.mapped_alleles_oppStrand # the complemnt of the new alleles # C  T
paste "_tmp.mapped_alleles" "_tmp.mapped_alleles_oppStrand" | \
    awk 'BEGIN{FS=OFS="\t"}{if($4=="+")print $1,$2,$3,$5,$6;else print $1,$2,$3,$7,$8;}' > _tmp.mapped_alleles_final # the new alleles on the + strand  # 11  21962991  21962991_Curly_f_ilmndup1   G   A  
paste _tmp.mapped_alleles_final  _tmp.ref_alleles  | awk 'BEGIN{FS=OFS="\t"}{if($4!=$6)print}' | wc -l ## 251 (all mismatches including indels) -- 173 if we used the known EquCab3 markers only
paste _tmp.mapped_alleles_final  _tmp.ref_alleles  | awk 'BEGIN{FS=OFS="\t"}{if($4!=$6)print}' | awk 'BEGIN{FS=OFS="\t"}{if($4=="-" || $5=="-")next;print}' | wc -l ## 134 (SNP mistmaches) -- 78 if we used the known EquCab3 markers only
grep "^#" _ref.vcf > _matchingSNPs.vcf
paste _tmp.mapped_alleles_final  _tmp.ref_alleles  | awk 'BEGIN{FS=OFS="\t"}{if($4==$6 && $5!="-")print $1,$2,$3,$4,$5,".",".","."}' >>  _matchingSNPs.vcf ##  $5!="-" remove additional 22 indels 
grep -v "^#" _matchingSNPs.vcf | wc -l ## 81672
grep -v "^#" _matchingSNPs.vcf | cut -f1 | sort | uniq -c | grep Un_ | wc -l ## 349



#snp="AX-102955572"
#paste _tmp.mapped_alleles_final  _tmp.ref_alleles  | awk 'BEGIN{FS=OFS="\t"}{if($2!=$4)print}' | grep -A10 $snp
#grep $snp Equine80select_24_20067593_B1_remapped.csv
#grep $snp temp_topseq.sam
#grep $snp temp_probe.sam

##################################
## Assessment of markers with the same position but different probe sequences or ref/alt alleles (for the REF matching markers only)
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11,$6}' | sort | uniq | wc -l ## $10 (chr),$11 (MapInfo/position),$6 (AlleleA_ProbeSeq) # 76235
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11,$23,$24}' | sort | uniq | wc -l ## $10 (chr),$11 (MapInfo/position), $23 (Ref_EquCab3),$24 (Alt_EquCab3) #76171
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11}' | sort | uniq | wc -l ## 76027

## same positions with different probes
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11,$6}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2}' | sort | uniq -c | awk '{if($1>1)print}' | sort -k1,1nr ## 69 positions with different probes

## same positions with different ref/alt alleles
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11,$23,$24}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2}' | sort | uniq -c | awk '{if($1>1)print}' | sort -k1,1nr ## 5 positions with different ref/alt alleles (same ref allele but different alternative)
#      2 1,109211964
#      2 16,21551060
#      2 25,27219807
#      2 3,79544174
#      2 3,79579925

## same positions with same probes but different ref/alt alleles
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11,$6,$23,$24}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2,$3}' | sort | uniq -c | awk '{if($1>1)print}' | sort -k1,1nr ## 2 positions with different ref/alt alleles for the same probe
#      2 25,27219807,TCATCGTCTTCTGGAGGAGAAGGTATCATGGAACTCTGAGATCCAGACTG
#      2 3,79544174,GGCTTTCTTTTCTCCCCCTCTCTCCTAATAGTGTATTCATAGGGACTTGG
grep '25,27219807' Equine80select_24_20067593_B1_remapped.csv | grep 'TCATCGTCTTCTGGAGGAGAAGGTATCATGGAACTCTGAGATCCAGACTG' ## Two different Alt alleles for the same reference allele
grep '3,79544174' Equine80select_24_20067593_B1_remapped.csv | grep 'GGCTTTCTTTTCTCCCCCTCTCTCCTAATAGTGTATTCATAGGGACTTGG' ## Two different Alt alleles for the same reference allele

## Remove likely polymorphic probes 
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=OFS=","}{print $10,$11,$23,$24}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2}' | sort | uniq -c | awk '{if($1>1)print $2}' > polymorphic_positions.txt
grep -Ff polymorphic_positions.txt Equine80select_24_20067593_B1_remapped.csv | cut -d, -f2 | grep -vFwf - _matchingSNPs.vcf > _matchingSNPs_binary.vcf
grep -v "^#" _matchingSNPs_binary.vcf | wc -l ## 81660
grep -v "^#" _matchingSNPs_binary.vcf | cut -f1 | sort | uniq -c | grep Un_ | wc -l ## 349
##################################
## Remove low quality probes markers 
cat inconsistent_probes | grep -vFwf - _matchingSNPs_binary.vcf > _matchingSNPs_binary_consistantMapping.vcf
grep -v "^#" _matchingSNPs_binary_consistantMapping.vcf | wc -l ## 80197
grep -v "^#" _matchingSNPs_binary_consistantMapping.vcf | cut -f1 | sort | uniq -c | grep Un_ | wc -l ## 116
grep -v "^#" _matchingSNPs_binary_consistantMapping.vcf | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,"0",$2,$4,$5}' | sort -k1,1d -k4,4n > Equine80select_remapped_equCab3.bim ## 80223
######################################################
## Create final map 
## Note: No need for EquCab2 coordinates because we can map byt SNP_ids (For EquCab2 coordinates, merge with the manifest)
## chr \t pos \t snpID \t SNP_alleles \t genomic_alleles \t SNP_ref_alleles \t genomic_ref_allele
## This map allows remapping to the same SNP alleles or PosStrand_alleles (i.e., VCF alleles). 
## In either case, their is a ref_allele to use in PLINK2
grep -v "^#" _matchingSNPs_binary_consistantMapping.vcf | cut -f3 > _tmp.matchingSNPs_binary_consistantMapping.snpIDs
tail -n+2 "$manifest" | cut -d, -f2,4 | awk 'BEGIN{FS=",";OFS="\t"}{gsub(/\[|\]/, "", $2);gsub(/\//, OFS, $2);print $0}' > _tmp.snp_alleles # the SNP name and alleles # 21962991_Curly_f_ilmndup1   A   G
cat _tmp.mapped_alleles_final | awk 'BEGIN{FS=OFS="\t"}{print $4,$5}' | tr 'TCGA' 'AGCT' > _tmp.mapped_alleles_final_complementary
paste _tmp.snp_alleles _tmp.mapped_alleles_final _tmp.mapped_alleles_final_complementary | grep -Fwf _tmp.matchingSNPs_binary_consistantMapping.snpIDs | awk 'BEGIN{FS=OFS="\t"}{\
    if($1==$6 && $2==$7 && $3==$8) print $4,$5,$6,$2","$3,$7","$8,$7,$7;  
    else if($1==$6 && $2==$8 && $3==$7) print $4,$5,$6,$2","$3,$8","$7,$7,$7;  
    else if($1==$6 && $2==$9 && $3==$10) print $4,$5,$6,$2","$3,$7","$8,$9,$7;  
    else if($1==$6 && $2==$10 && $3==$9) print $4,$5,$6,$2","$3,$8","$7,$9,$7;  
    else print "Error",$0}' > matchingSNPs_binary_consistantMapping.map
grep "^Error" matchingSNPs_binary_consistantMapping.map | wc -l ## 0


#################################

# Convert Equine80select manifest to Plink BIM file for equCab2 and equCab3 separately with Allele 1 and Allele 2 correspond to REF and ALT alleles respectively.
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=",";OFS="\t"}{if($10 && $9=="2")print $10,$2,"0",$11}' | sort -k2,2 > Equine80select2.map  ## 75952
tail -n+2 Equine80select_24_20067593_B1_remapped.csv | awk 'BEGIN{FS=",";OFS="\t"}{if($10 && $9=="3")print $10,$2,"0",$11}' | sort -k2,2 > Equine80select3.map ## 6022



#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc670k.unique_remap.FINAL.csv.gz
#gunzip MNEc670k.unique_remap.FINAL.csv.gz
mollyFile=$parentageDir/backup_original/MNEc670k.unique_remap.FINAL.csv
cat $mollyFile | awk -F"," '{if($10)print}'  | sed 's/chrUn_ref|/Un_/' | sed 's/\.1|,/v1,/' | sed 's/,chr/,/g' | tr ',' '\t' > MNEc670k_remap.tab
mollyMap="$work_dir"/MNEc670k_remap.tab  ## 641919
awk 'BEGIN{FS=OFS="\t"}{print $5"."$7,$6"."$8,$10}' $mollyMap | head -n2
#EC2_chrom.EC2_pos       EC3_chrom.EC3_pos       EC3_ALT
#1.3745                  1.3052                  C
awk 'BEGIN{FS=OFS="\t"}{print $5"."$7}' $mollyMap | sort | uniq -c | awk '{if($1>1)print}' ## no duplicates
awk 'BEGIN{FS=OFS="\t"}{print $6"."$8}' $mollyMap | sort | uniq -c | awk '{if($1>1)print}' ## no duplicates

for f in Equine80select2.map Equine80select3.map;do
  f2=${f%.map}_molly_equCab3.bim
  equ2=$(comm -12 <(cat $f | awk '{print $1"."$4}' | sort) <(cat $mollyMap |  awk '{print $5"."$7}' | sort) | wc -l)
  equ3=$(comm -12 <(cat $f | awk '{print $1"."$4}' | sort) <(cat $mollyMap |  awk '{print $6"."$8}' | sort) | wc -l)
  if [ $equ2 -gt $equ3 ];then
    echo "equ2=$equ2 & equ3=$equ3 ... $f is equ2"
    awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($10){a[$5"_"$7]=$6;b[$5"_"$7]=$8;c[$5"_"$7]=$9 FS $10;}next}{if(a[$1"_"$4])print a[$1"_"$4] FS $2 FS $3 FS b[$1"_"$4] FS c[$1"_"$4];}' "$mollyMap" $f | sort -k1,1d -k4,4n > $f2
  else
    echo "equ2=$equ2 & equ3=$equ3 ... $f is equ3"
    awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($10){a[$6"_"$8]=$6;b[$6"_"$8]=$8;c[$6"_"$8]=$9 FS $10;}next}{if(a[$1"_"$4])print a[$1"_"$4] FS $2 FS $3 FS b[$1"_"$4] FS c[$1"_"$4];}' "$mollyMap" $f | sort -k1,1d -k4,4n > $f2
  fi
done
# equ2=66209 & equ3=26 ... Equine80select2.map is equ2
# equ2=3 & equ3=4955 ... Equine80select3.map is equ3

## merge the bim files of Equine80select
cat Equine80select3_molly_equCab3.bim Equine80select2_molly_equCab3.bim | sort -k1,1d -k4,4n > Equine80select_molly_equCab3.bim 
wc -l *_molly_equCab3.bim 
#  70987 Equine80select2_molly_equCab3.bim
#   5162 Equine80select3_molly_equCab3.bim
#  76149 Equine80select_molly_equCab3.bim

# Compare the two remapping results
comm -12 <(cat Equine80select_molly_equCab3.bim | cut -d$'\t' -f2 | sort) <(cat Equine80select_remapped_equCab3.bim | cut -d$'\t' -f2 | sort) | wc -l  ## 75684 shared 
comm -12 <(cat Equine80select_molly_equCab3.bim | cut -d$'\t' -f2 | sort) <(cat Equine80select_remapped_equCab3.bim | cut -d$'\t' -f2 | sort) > molly_remapped_common_snps.txt
comm -12 <(cat molly_remapped_common_snps.txt | grep -Fwf - Equine80select_molly_equCab3.bim | sort) <(cat molly_remapped_common_snps.txt | grep -Fwf - Equine80select_remapped_equCab3.bim | sort) | wc -l  ## 76053
comm -12 <(cat molly_remapped_common_snps.txt | grep -Fwf - Equine80select_molly_equCab3.bim | sort) <(cat molly_remapped_common_snps.txt | grep -Fwf - Equine80select_remapped_equCab3.bim | sort) | cut -d$'\t' -f2 > molly_remapped_matching_snps.txt
paste <(cat molly_remapped_matching_snps.txt | grep -vFwf - molly_remapped_common_snps.txt | grep -Fwf - Equine80select_molly_equCab3.bim | sort -k2,2) <(cat molly_remapped_matching_snps.txt | grep -vFwf - molly_remapped_common_snps.txt | grep -Fwf - Equine80select_remapped_equCab3.bim | sort -k2,2) | head

"""
Un_NW_019642777v1       BIEC2_1009547           0       4472            C       A       Un_NW_019642671v1       BIEC2_1009547           0       4871            A       C   swabbed alleles in Molly files
11                      BIEC2_156113            0       30686647        A       C       11                      BIEC2_156113            0       30648596        C       A   swabbed alleles in Molly file
13                      BIEC2_204488            0       3528848         T       C       13                      BIEC2_204488            0       3528858         T       C
13                      BIEC2_226662            0       26373342        C       T       13                      BIEC2_226662            0       26521282        C       T
14                      BIEC2_247157            0       16998719        C       T       8                       BIEC2_247157            0       40773756        T       C   swabbed alleles in Molly file
15                      BIEC2_289793_ilmndup1   0       15955711        C       T       15                      BIEC2_289793_ilmndup1   0       16032903        C       T
15                      BIEC2_289793_ilmndup2   0       15955711        C       T       15                      BIEC2_289793_ilmndup2   0       15994103        C       T
16                      BIEC2_356910            0       37385819        C       A       16                      BIEC2_356910            0       37347586        A       C   swabbed alleles in Molly file
20                      BIEC2_565175            0       51490766        C       T       20                      BIEC2_565175            0       51352858        C       T   
26                      BIEC2_696038            0       39707112        G       T       26                      BIEC2_696038            0       39784486        A       C   swabbed and flipped alleles in Molly file
"""
snp="Affx-101373202"
#paste _tmp.mapped_alleles_final  _tmp.ref_alleles  | awk 'BEGIN{FS=OFS="\t"}{if($2!=$4)print}' | grep -A10 $snp
grep $snp Equine80select_24_20067593_B1_remapped.csv
grep $snp temp_topseq.sam
grep $snp temp_probe.sam

