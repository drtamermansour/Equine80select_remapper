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


