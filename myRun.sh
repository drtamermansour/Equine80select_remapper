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


##############################################################################################################################
## This section was intended to test the output of remap_manifest.py to assess its performance and decide the best filters 
#############################################################################################################################
## Explore the manifest header and format 
mkdir -p temp_explore
outManifest="$work_dir"/results/Equine80select_24_20067593_B1_remapped_equCab3.csv
grep IlmnID $origManifest | tr ',' '\n' > temp_explore/1.txt
grep IlmnID $outManifest | tr ',' '\n' | awk '{print NR}' > temp_explore/2.txt
grep IlmnID $outManifest | tr ',' '\n' > temp_explore/3.txt
grep -A1 IlmnID $outManifest | tail -n1 | tr ',' '\n' > temp_explore/4.txt
paste temp_explore/1.txt temp_explore/2.txt temp_explore/3.txt temp_explore/4.txt

## Find the possiblities of output Strand_EquCab3 ($22)
tail -n+2 $outManifest | awk -F, '{print $22}' | sort | uniq -c
##  40521 +
##  41424 -
##     29 N/A ## These records has N in the output Ref_EquCab3 & Alt_EquCab3 ($23 & $24); likely failed to align 

##########################
## Generate histogram of MAPQ scores for all TopGenomicSeq alignments
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($25/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $outManifest) > temp_explore/MAPQ_TopGenomicSeq.histo

## Generate histogram of MAPQ scores for all probes alignments
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($26/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $outManifest) > temp_explore/MAPQ_Probe.histo

cat <(grep -v "^@" results/temp/temp_topseq.sam | sed 's/_A\t/\t/;s/_B\t/\t/;') <(grep -v "^@" results/temp/temp_probe.sam) | cut -f1,3 | sort | uniq -c | awk '{if($1!=3)print}' > temp_explore/inconsistant_alignments 
cat temp_explore/inconsistant_alignments | awk '{print $2}' | sort | uniq > temp_explore/inconsistent_probes
cat temp_explore/inconsistent_probes | grep -Fwf - $outManifest > temp_explore/inconsistent_remapped.csv

## Generate histogram of MAPQ scores for TopGenomicSeq alignments from inconsistent_probes
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($25/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 temp_explore/inconsistent_remapped.csv) > temp_explore/MAPQ_TopGenomicSeq_inconsistent_remapped.histo

## Generate histogram of MAPQ scores for probes with inconsistent alignments
awk -v size=2 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($26/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 temp_explore/inconsistent_remapped.csv) > temp_explore/MAPQ_Probe_inconsistent_remapped.histo


##################################
## Assessment of markers with the same position but different probe sequences or ref/alt alleles (for the REF matching markers only)
tail -n+2 $outManifest | awk 'BEGIN{FS=OFS=","}{print $10,$11,$6}' | sort | uniq | wc -l ## $10 (chr),$11 (MapInfo/position),$6 (AlleleA_ProbeSeq) # 76235
tail -n+2 $outManifest | awk 'BEGIN{FS=OFS=","}{print $10,$11,$23,$24}' | sort | uniq | wc -l ## $10 (chr),$11 (MapInfo/position), $23 (Ref_EquCab3),$24 (Alt_EquCab3) #76172
tail -n+2 $outManifest | awk 'BEGIN{FS=OFS=","}{print $10,$11}' | sort | uniq | wc -l ## 76166

## same positions with different probes
tail -n+2 $outManifest | awk 'BEGIN{FS=OFS=","}{print $10,$11,$6}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2}' | sort | uniq -c | awk '{if($1>1)print}' | sort -k1,1nr ## 69 positions with different probes

## same positions with different ref/alt alleles
tail -n+2 $outManifest | awk 'BEGIN{FS=OFS=","}{print $10,$11,$23,$24}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2}' | sort | uniq -c | awk '{if($1>1)print}' | sort -k1,1nr ## 6 positions with different ref/alt alleles (same ref allele but different alternative)
#      2 1,109211964
#      2 16,21551060
#      2 25,27219807
#      2 3,79544174
#      2 3,79579925
#      2 6,52969965

## same positions with same probes but different ref/alt alleles
tail -n+2 $outManifest | awk 'BEGIN{FS=OFS=","}{print $10,$11,$6,$23,$24}' | sort | uniq | \
    awk 'BEGIN{FS=OFS=","}{print $1,$2,$3}' | sort | uniq -c | awk '{if($1>1)print}' | sort -k1,1nr ## 3 positions with different ref/alt alleles for the same probe
#      2 25,27219807,TCATCGTCTTCTGGAGGAGAAGGTATCATGGAACTCTGAGATCCAGACTG
#      2 3,79544174,GGCTTTCTTTTCTCCCCCTCTCTCCTAATAGTGTATTCATAGGGACTTGG
#      2 6,52969965,AACAGCTGATTATGTTTCAGCGGGACCATCTGTTCGAAACCCCAAAGCAC
grep '25,27219807' $outManifest  | grep 'TCATCGTCTTCTGGAGGAGAAGGTATCATGGAACTCTGAGATCCAGACTG' ## Two different Alt alleles for the same reference allele
grep '3,79544174' $outManifest  | grep 'GGCTTTCTTTTCTCCCCCTCTCTCCTAATAGTGTATTCATAGGGACTTGG' ## Two different Alt alleles for the same reference allele


#############################################################################
## Prep. work to decide the best way to create "allele_usage_decision.txt"
#############################################################################
## annotate records where SNP ($4) != alleles reported in TopGenomicSeq ($18)
cat $outManifest | sed 's/[^,]*\(\[[^]]*\]\)[^,]*/\1/g' > temp_explore/annotated_discrepancies_manifest.csv
cat temp_explore/annotated_discrepancies_manifest.csv | awk 'BEGIN{FS=OFS=","}{print $3,$16,$22,$4,$17,$18}' > temp_explore/check_output.csv
tail -n+2 temp_explore/check_output.csv | sort | uniq -c > temp_explore/check_output.groups
## Example output
##    269 BOT,BOT,+,[C/T],[A/G],[A/G]
##      3 BOT,BOT,N/A,[G/T],[A/C],[A/C]
##   4360 BOT,BOT,-,[G/T],[A/C],[A/C]
##  14510 BOT,TOP,+,[C/T],[A/G],[A/G]
##   3696 BOT,TOP,+,[G/T],[A/C],[A/C]
##      6 MINUS,PLUS,+,[D/I],[-/CAGAAAAGAAG],[A/C/G/T]


## Conclusion:
## The remapping script should do the following:
## 1. Exclude records with N/A in Strand_EquCab3 ($22) from further analysis
## 2. For now, we can discard records with indels (I/D) in SNPs ($4) as they are few and complex to handle
## 3. if the PLINK output is not using the same manifest alleles ($4), exclude ambiguous SNPs, then:
##    a. if SNP alleles ($4) match TopGenomicSeq alleles ($18) && Strand_EquCab3 ($22) == +, use SNP alleles as is
##    b. if SNP alleles ($4) match TopGenomicSeq alleles ($18) && Strand_EquCab3 ($22) == -, complement SNP alleles
##    c. if SNP alleles ($4) do not match TopGenomicSeq alleles ($18) && Strand_EquCab3 ($22) == +, complement SNP alleles
##    d. if SNP alleles ($4) do not match TopGenomicSeq alleles ($18) && Strand_EquCab3 ($22) == -, use SNP alleles as is
## Otherwise decide if the SNP alleles ($4) should be used as is or be complemented (without matching the alleles to avoid mistakes of ambiguos SNPs):
##    a. IlmnStrand ($3) == SourceStrand ($16) && SourceSeq ($17) == TopGenomicSeq ($18) && Strand_EquCab3 ($22) == +  => use SNP alleles as is
##    b. IlmnStrand ($3) == SourceStrand ($16) && SourceSeq ($17) == TopGenomicSeq ($18) && Strand_EquCab3 ($22) == -  => complement SNP alleles
##    c. IlmnStrand ($3) != SourceStrand ($16) && SourceSeq ($17) == TopGenomicSeq ($18) && Strand_EquCab3 ($22) == +  => complement SNP alleles
##    d. IlmnStrand ($3) != SourceStrand ($16) && SourceSeq ($17) == TopGenomicSeq ($18) && Strand_EquCab3 ($22) == -  => use SNP alleles as is
##    e. IlmnStrand ($3) == SourceStrand ($16) && SourceSeq ($17) != TopGenomicSeq ($18) && Strand_EquCab3 ($22) == +  => complement SNP alleles
##    f. IlmnStrand ($3) == SourceStrand ($16) && SourceSeq ($17) != TopGenomicSeq ($18) && Strand_EquCab3 ($22) == -  => use SNP alleles as is
##    g. IlmnStrand ($3) != SourceStrand ($16) && SourceSeq ($17) != TopGenomicSeq ($18) && Strand_EquCab3 ($22) == +  => use SNP alleles as is
##    h. IlmnStrand ($3) != SourceStrand ($16) && SourceSeq ($17) != TopGenomicSeq ($18) && Strand_EquCab3 ($22) == -  => complement SNP alleles

## Apply the above logic to summarize the number of records that should use SNP alleles as is or be complemented
tail -n+2 temp_explore/annotated_discrepancies_manifest.csv | awk -F, 'BEGIN{OFS="\t";a["as_is"]=0; a["complement"]=0; b=0} {
if($22=="N/A") next;
if($4=="[D/I]" || $4=="[I/D]") next;
if($3==$16 && $17==$18 && $22=="+") a["as_is"]+=1;
    else if($3==$16 && $17==$18 && $22=="-") a["complement"]+=1;
    else if($3!=$16 && $17==$18 && $22=="+") a["complement"]+=1;
    else if($3!=$16 && $17==$18 && $22=="-") a["as_is"]+=1;
    else if($3==$16 && $17!=$18 && $22=="+") a["complement"]+=1;
    else if($3==$16 && $17!=$18 && $22=="-") a["as_is"]+=1;
    else if($3!=$16 && $17!=$18 && $22=="+") a["as_is"]+=1;
    else if($3!=$16 && $17!=$18 && $22=="-") a["complement"]+=1;

if($3==$16 && $17==$18 && $28!=$29) {b+=1;print $0",st1"}
    else if($3!=$16 && $17==$18 && $30!=$29) {b+=1;print $0",st2"}
    else if($3==$16 && $17!=$18 && $30!=$29) {b+=1;print $0",st3"}
    else if($3!=$16 && $17!=$18 && $28!=$29) {b+=1;print $0",st4"}

} END {for(i in a) print i,a[i]; print "Total mismatches:",b}' 
## complement      37786
## as_is   44020
## Total mismatches:       0

tail -n+2 temp_explore/annotated_discrepancies_manifest.csv | awk -F, 'BEGIN{OFS="\t";a["as_is"]=0; a["complement"]=0; b=0} {
if($22=="N/A") next;
if(($4=="[D/I]" || $4=="[I/D]") && $3==$16 && $17==$18 && $22=="+") print $2,"indel_as_is";
else if(($4=="[D/I]" || $4=="[I/D]") && $3==$16 && $17==$18 && $22=="-") print $2,"indel_complement";
else if(($4=="[D/I]" || $4=="[I/D]") && $3!=$16 && $17==$18 && $22=="+") print $2,"indel_complement";
else if(($4=="[D/I]" || $4=="[I/D]") && $3!=$16 && $17==$18 && $22=="-") print $2,"indel_as_is";
else if(($4=="[D/I]" || $4=="[I/D]") && $3==$16 && $17!=$18 && $22=="+") print $2,"indel_complement";
else if(($4=="[D/I]" || $4=="[I/D]") && $3==$16 && $17!=$18 && $22=="-") print $2,"indel_as_is";
else if(($4=="[D/I]" || $4=="[I/D]") && $3!=$16 && $17!=$18 && $22=="+") print $2,"indel_as_is";
else if(($4=="[D/I]" || $4=="[I/D]") && $3!=$16 && $17!=$18 && $22=="-") print $2,"indel_complement";
else if($3==$16 && $17==$18 && $22=="+") print $2,"as_is";
else if($3==$16 && $17==$18 && $22=="-") print $2,"complement";
else if($3!=$16 && $17==$18 && $22=="+") print $2,"complement";
else if($3!=$16 && $17==$18 && $22=="-") print $2,"as_is";
else if($3==$16 && $17!=$18 && $22=="+") print $2,"complement";
else if($3==$16 && $17!=$18 && $22=="-") print $2,"as_is";
else if($3!=$16 && $17!=$18 && $22=="+") print $2,"as_is";
else if($3!=$16 && $17!=$18 && $22=="-") print $2,"complement";
}' > temp_explore/allele_usage_decision.txt


###############################################################
## Benchmark aganist the already known EquCab3 markers 
###############################################################
mkdir -p remap_assessment
## Summary statistics of input EquCab2 and EquCab3 markers of the original manifest
cat $outManifest | awk 'BEGIN{FS=OFS=","}{print $9}' | sort | uniq -c
#  75952 2.0
#   6022 3.0

## Summary statistics of remapping results ($20 is the new Chr_EquCab3)
cat $outManifest | awk 'BEGIN{FS=OFS=","}{print $20}' | sort | uniq -c | awk '{if($2==0) print}'  ## 29 with 0 value (i.e. failed to remap)
cat $outManifest | awk 'BEGIN{FS=OFS=","}{print $20}' | grep Un | wc -l ## 648 assigned to chrUn


#$9=Original_Assembly (2=EquCab2, 3=EquCab3)
#$2=Name
#$4=SNP e.g., [A/G] 
#$10:$11=Chr:MapInfo
#$3/$16=IlmnStrand/SourceStrand
#$20:$21=Chr_EquCab3:MapInfo_EquCab3
#$22=Strand_EquCab3

## confirm that input Chr:MapInfo ($10:$11) matched the output Chr_EquCab3:MapInfo_EquCab3 ($20:$21) for the known EquCab3 markers only ($9==3)
cat $outManifest | awk 'BEGIN{FS=OFS=","}{if($9==3)print $2,$4,$10":"$11,$3"/"$16,$20":"$21,$22}' | awk 'BEGIN{FS=","}{if($3==$5)print}' | wc -l ## 5934 remapped correctly
cat $outManifest | awk 'BEGIN{FS=OFS=","}{if($9==3)print $2,$4,$10":"$11,$3"/"$16,$20":"$21,$22}' | awk 'BEGIN{FS=","}{if($3!=$5)print}' > remap_assessment/equCab3.mismatches
cat remap_assessment/equCab3.mismatches | wc -l ## 88
cat remap_assessment/equCab3.mismatches | awk 'BEGIN{FS=OFS=","}{print $2,$4,$6}' | sort | uniq -c
#      9 [D/I],MINUS/PLUS,+
#      6 [I/D],MINUS/PLUS,+
#     15 [I/D],PLUS/PLUS,+
#      1 [A/C],TOP/TOP,+
#      5 [A/G],TOP/BOT,+
#      7 [A/G],TOP/BOT,-
#      4 [A/G],TOP/TOP,+
#      3 [A/G],TOP/TOP,-
#      8 [A/T],TOP/BOT,-
#      1 [C/G],TOP/TOP,+
#      1 [C/G],TOP/TOP,-
#      3 [G/C],BOT/TOP,+
#      2 [G/C],BOT/TOP,-
#      4 [T/C],BOT/BOT,+
#     11 [T/C],BOT/BOT,-
#      3 [T/C],BOT/TOP,-
#      1 [T/G],BOT/BOT,+
#      1 [T/G],BOT/BOT,-
#      1 [T/G],BOT/TOP,+
#      2 [T/G],BOT/TOP,-


#####################################
## Molly Cross-Validation
#####################################
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
comm -12 <(cat molly_remapped_common_snps.txt | grep -Fwf - Equine80select_molly_equCab3.bim | sort) <(cat molly_remapped_common_snps.txt | grep -Fwf - Equine80select_remapped_equCab3.bim | sort) | wc -l  ## 75649 shared with same SNP ID, chr and pos
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

