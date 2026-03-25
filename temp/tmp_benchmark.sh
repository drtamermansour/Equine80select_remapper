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
