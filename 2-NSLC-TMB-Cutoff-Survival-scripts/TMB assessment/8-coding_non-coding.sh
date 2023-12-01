#!bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Annotation of coding and non-coding variants using ANNOVAR. 
# There are two options possible using output of script 5-final_filtering.sh : 1) annotation of variants gathered from all three variant calling tools or 2) annotation of variants common to at least two variant calling tools.
# The annotation results are saved in a csv file.

Patient=$1

if [ "$Patient" = "" ]; then
echo "Error: missing patient identifier argument."
exit
fi

cd /mnt/raid5_hdd/ruelou01/NSLC_${Patient}

# Option 1 ******************************************************************************************************************
# Coding and non-coding annotation of filtered common SNVs and indels **to all variant calling tools used**

#/mnt/raid5_hdd/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo common_snv_indel.vcf --outfile common_snv_indel.avinput

#/mnt/raid5_hdd/annovar/annotate_variation.pl -build hg19 -out annotation/coding_non-coding_${Patient} common_snv_indel.avinput /mnt/raid5_hdd/annovar/humandb/

#rm common_snv_indel.avinput


# Option 2 ******************************************************************************************************************
# Coding and non-coding annotation of filtered common SNVs and indels **to at least two variant calling tools used**
# Removing unmapped contigs variants starting with GL000
grep -v "GL000" NSLC_${Patient}_common_2_snv_indel.vcf > NSLC_${Patient}_common_2_snv_indel_temp.vcf
cat NSLC_${Patient}_common_2_snv_indel_temp.vcf > NSLC_${Patient}_common_2_snv_indel.vcf

/mnt/raid5_hdd/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo NSLC_${Patient}_common_2_snv_indel.vcf --outfile common_2_snv_indel.avinput

/mnt/raid5_hdd/annovar/annotate_variation.pl -build hg19 -out annotation/coding_non-coding_${Patient} NSLC_${Patient}_common_2_snv_indel.avinput /mnt/raid5_hdd/annovar/humandb/

rm NSLC_${Patient}_common_2_snv_indel_temp.vcf
rm NSLC_${Patient}_common_2_snv_indel.avinput


#******************************************************************************************************************
# Table creation of the resulting annotation for specific patients
printf '"chr";"start";"end";"REF";"ALT";"Func.refGene";"Gene";"Mutation_type"' > annotation/coding_non-coding_${Patient}.csv
printf "\n" >> annotation/coding_non-coding_${Patient}.csv

# This awk function reads the exonic annotation file (*.exonic_variant_function) first followed by the full annotation file (*.variant_function) and extract the gene and mutation type of each variant. Afterwards, it registers non-coding or coding annotation of that variant in a csv file. If the variant is a coding variant, it also registers the corresponding full gene name and mutation type found in the exonic annotation file.
awk -F"\t" 'BEGIN{OFS="\";\""}; NR==FNR{gene[$14]=$3;mutation[$14]=$2;next;} NR!=FNR{if($1!="exonic"){print "\""$3,$4,$5,$6,$7,$1,$2"\""}; if($1=="exonic"){print "\""$3,$4,$5,$6,$7,$1,gene[$13],mutation[$13]"\""}}' annotation/coding_non-coding_${Patient}.exonic_variant_function annotation/coding_non-coding_${Patient}.variant_function >> annotation/coding_non-coding_${Patient}.csv



