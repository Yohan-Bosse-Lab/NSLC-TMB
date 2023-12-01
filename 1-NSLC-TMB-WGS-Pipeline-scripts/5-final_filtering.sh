#!/bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Sorting variants that are common to multiple tools to reduce false positives variants found in the vcfs of variant calling tools.

# NOTE: There is potentially a lot of unecessary files (such as extra vcf files) that are not used inside this script or by any folowing scripts. Some optimization is needed.

Patient=$1

if [ "$Patient" = "" ]; then
echo "Error: missing patient identifier argument."
exit
fi

cd /mnt/raid5_hdd/ruelou01/NSLC_${Patient}

# Removing Strelka file containing pre database filtered variants without the added GT tag. The GT tag is required for the later annotation in script 8-coding_non-coding.sh". In short, we keep "Strelka/${Patient}_Strelka_preFilter_GT.vcf" and we remove "Strelka/${Patient}_Strelka_preFilter.vcf".
rm Strelka/${Patient}_Strelka_preFilter.vcf

## This is the database filtering: extracting variants with MAF (minor allele frequency) < 0.001 (following the method from Zhang et al. 2021) from databases for both formal and inverted REF and ALT columns order.
## Formal column order:
# The first command extract only the variant data lines. The second command in the pipeline extracts the ID of the variant (column 3) and the database annotation information (column 8). The third command removes all character contained in the column 8 and extracts the MAF values of 1000 Genome, ExAC and gnomAD databases, in that order. The fourth command 1) replaces "." entries in MAF fields as 0; 2) filters the MAF values of variants and keeps those < 0.001 for all three databases. Finally, it is outputed as a temporary file to extract the full records of the remaining (filtered) variants in the next step.
grep -v '#' annotation/formal_annotation.hg19_multianno.vcf | awk '{print $3, $8}' | sed -e 's/ .*1000g2015aug_all=/ /g' -e 's/;ExAC_ALL=/ /g' -e 's/;.*AF=/ /g' -e 's/;.*//g' | awk '{for(i=1; i<=NF; i++) {if($i==".") $i=0}} {if($2<0.001 && $3<0.001 && $4<0.001) print}' > annotation/formal_db_filtering_results.txt

## These steps are for extracting full records of the database filtered variants and make the new vcf file containing only the database filtered variants.
# Extracting the databse filtered variants ID.
awk '{print $1}' annotation/formal_db_filtering_results.txt > annotation/formal_ID.temp
# Extracting all corresponding filtered variants from the full records.
grep '#' annotation/formal_annotation.hg19_multianno.vcf > annotation/formal_dbFilter.vcf
grep -Fwf annotation/formal_ID.temp annotation/formal_annotation.hg19_multianno.vcf >> annotation/formal_dbFilter.vcf


## Inverted column order:
# Same principles as the formal column order (see above).
grep -v '#' annotation/inverted_annotation.hg19_multianno.vcf | awk '{print $3, $8}' | sed -e 's/ .*1000g2015aug_all=/ /g' -e 's/;ExAC_ALL=/ /g' -e 's/;.*AF=/ /g' -e 's/;.*//g' | awk '{for(i=1; i<=NF; i++) {if($i==".") $i=0}} {if($2<0.001 && $3<0.001 && $4<0.001) print}' > annotation/inverted_db_filtering_results.txt

## These steps are for extracting full records of the database filtered variants and make the new vcf file containing only the database filtered variants.
# Extracting the databse filtered variants ID.
awk '{print $1}' annotation/inverted_db_filtering_results.txt > annotation/inverted_ID.temp
# Extracting all corresponding filtered variants from the full records.
grep '#' annotation/inverted_annotation.hg19_multianno.vcf > annotation/inverted_dbFilter.vcf
grep -Fwf annotation/inverted_ID.temp annotation/inverted_annotation.hg19_multianno.vcf >> annotation/inverted_dbFilter.vcf


for tool in TNsnv Mutect2 TNscope Strelka
do
# Extracting tool variants from database filtering into a simple ID text file (for Venn diagram purposes) and a vcf file (for further analysis). The extraction is done for both formal and inverted REF and ALT columns order.

# Formal column order:
# ID text file (for Venn diagram) and excluding unmapped variants (from contigs starting with 'GL000'):
grep -v '#' annotation/formal_dbFilter.vcf | grep -v 'GL000' | awk '{print $3}' | grep "${tool}" > ${tool}/${Patient}_${tool}_formal_dbFilter_ID.txt
# vcf file (for further analysis) and excluding unmapped variants (from contigs starting with 'GL000'):
grep '#' annotation/formal_dbFilter.vcf > ${tool}/${Patient}_${tool}_formal_dbFilter.vcf
grep -v '#' annotation/formal_dbFilter.vcf | grep -v 'GL000' | grep "${tool}" >> ${tool}/${Patient}_${tool}_formal_dbFilter.vcf

# Inverted column order:
# ID text file (for Venn diagram) and excluding unmapped variants (from contigs starting with 'GL000'):
grep -v '#' annotation/inverted_dbFilter.vcf | grep -v 'GL000' | awk '{print $3}' | grep "${tool}" > ${tool}/${Patient}_${tool}_inverted_dbFilter_ID.txt
# vcf file (for further analysis) and excluding unmapped variants (from contigs starting with 'GL000'):
grep '#' annotation/inverted_dbFilter.vcf > ${tool}/${Patient}_${tool}_inverted_dbFilter.vcf
grep -v '#' annotation/inverted_dbFilter.vcf | grep -v 'GL000' | grep "${tool}" >> ${tool}/${Patient}_${tool}_inverted_dbFilter.vcf

# Quick sorting the ID files (for Venn diagram).
sort ${tool}/${Patient}_${tool}_formal_dbFilter_ID.txt > ${tool}/${Patient}_${tool}_f_temp.txt
sort ${tool}/${Patient}_${tool}_inverted_dbFilter_ID.txt > ${tool}/${Patient}_${tool}_i_temp.txt

# Extracting common variants that have MAF <0.001 (both in formal or inverted REF and ALT orders), which means variants that are both in the formal and inverted filtered files. The identification of such variants is done by using the ID column, which is the same for formal order or inverted order.
comm -12 ${tool}/${Patient}_${tool}_f_temp.txt ${tool}/${Patient}_${tool}_i_temp.txt | sort -V | awk -v Patient="$Patient" -v tool="$tool" 'BEGIN {OFS=FS = ":"} ; {if (length($4)==1 && length($5)==1){print $2, $3, $4, $5 > tool"/"Patient"_"tool"_completeFilter_SNP.txt"}; if (length($4)!=length($5)){print $2, $3, $4, $5 > tool"/"Patient"_"tool"_completeFilter_Ind.txt"}}'

#rm ${tool}/${Patient}_${tool}_f_temp.txt
#rm ${tool}/${Patient}_${tool}_i_temp.txt

done

#####################################################################################################

# Get common filtered SNPs for TNsnv, Mutect2 and TNscope
sort TNsnv/${Patient}_TNsnv_completeFilter_SNP.txt > sort_TNsnv.SNP
sort Mutect2/${Patient}_Mutect2_completeFilter_SNP.txt > sort_Mutect2.SNP
sort TNscope/${Patient}_TNscope_completeFilter_SNP.txt > sort_TNscope.SNP

comm -12 sort_TNsnv.SNP sort_Mutect2.SNP > SNPs/filtered_TNsnv_Mutect2_SNP.txt
comm -12 sort_TNsnv.SNP sort_TNscope.SNP > SNPs/filtered_TNsnv_TNscope_SNP.txt
comm -12 sort_Mutect2.SNP sort_TNscope.SNP > SNPs/filtered_Mutect2_TNscope_SNP.txt
comm -12 SNPs/filtered_TNsnv_Mutect2_SNP.txt sort_TNscope.SNP > common_SNV.txt

# Creating another file for common SNPs to at least 2 tools. See below for further details.
cat SNPs/filtered_TNsnv_Mutect2_SNP.txt SNPs/filtered_TNsnv_TNscope_SNP.txt SNPs/filtered_Mutect2_TNscope_SNP.txt | sort -V | uniq > common_2_SNV.txt

# Get common filtered indels for Mutect2, TNscope and Strelka
sort Mutect2/${Patient}_Mutect2_completeFilter_Ind.txt > sort_Mutect2.Ind
sort TNscope/${Patient}_TNscope_completeFilter_Ind.txt > sort_TNscope.Ind
sort Strelka/${Patient}_Strelka_completeFilter_Ind.txt > sort_Strelka.Ind

comm -12 sort_Mutect2.Ind sort_TNscope.Ind > Indels/filtered_Mutect2_TNscope_Ind.txt
comm -12 sort_Mutect2.Ind sort_Strelka.Ind > Indels/filtered_Mutect2_Strelka_Ind.txt
comm -12 sort_TNscope.Ind sort_Strelka.Ind > Indels/filtered_TNscope_Strelka_Ind.txt
comm -12 Indels/filtered_Mutect2_TNscope_Ind.txt sort_Strelka.Ind > common_Ind.txt

# Creating another file for SNPs common to at least 2 tools. See below for further details.
cat Indels/filtered_Mutect2_TNscope_Ind.txt Indels/filtered_Mutect2_Strelka_Ind.txt Indels/filtered_TNscope_Strelka_Ind.txt | sort -V | uniq > common_2_Ind.txt

#####################################################################################################

## ***IF ONLY VARIANTS COMMON TO ALL 3 TOOLS ARE TAKEN FOR TMB***

# Extract common SNPs and indels from the Mutect2 reference vcf and create a vcf ready for further annotation. See "8-coding_non-coding.sh" for further details.
# The choice of Mutect2 as the reference is arbitrary and has no impact on the created vcf (since, at this point, the output variants are the same across TNsnv, Mutect2, TNscope and Strelka).
#check=0
#for line in $(cat common_SNV.txt common_Ind.txt | sort -V)
#do
#if [ $check -eq 0 ]; then
#grep "#" Mutect2/${Patient}_Mutect2_preFilter.vcf > common_snv_indel.temp
#check=1
#fi
#grep "${line}" Mutect2/${Patient}_Mutect2_preFilter.vcf >> common_snv_indel.temp
#done

#sed -e 's/Mutect2_//g' -e 's/Mutect2://g' -e 's/Strelka://g' -e 's/TNscope://g' -e 's/TNsnv://g' common_snv_indel.temp > common_snv_indel.vcf

#####################################################################################################

## ***IF VARIANTS COMMON TO AT LEAST 2 TOOLS ARE TAKEN FOR TMB***

# Puts together filtered SNVs and indels.
cat common_2_SNV.txt common_2_Ind.txt | sort -V > snvs_ind.temp

# Here, we just take the vcf header on the patient's Mutect2 file. The header only include line starting with #. Need to do this so the final file (with common variants to 2 tools) has a valid vcf format. This could have been taken from any other tool (TNscope, TNsnv or Strelka).
grep "#" Mutect2/${Patient}_Mutect2_preFilter.vcf > NSLC_${Patient}_common_2_snv_indel.vcf

# The following pipeline extracts common SNPs and indels to at least 2 tools and create a vcf ready for further annotation. See "8-coding_non-coding.sh" for further details.
# The first command matches filtered variants with the previous complete ("complete" means not set in a formal or inverted order) vcf file specific their tools. The second command remove the tool tags found in the ID column, because we need to remove duplicates (e.g. if Mutect2 and TNscope have the same filtered variants, we don't want a duplicated record of that variant for the final analyses). The third command remove duplicated filtered variant record.
grep -h -Fwf snvs_ind.temp */${Patient}_*_preFilt*.vcf | sed -e 's/Mutect2_//g' -e 's/Mutect2://g' -e 's/Strelka://g' -e 's/TNscope://g' -e 's/TNsnv://g' | sort -u -V -k3,3 >> NSLC_${Patient}_common_2_snv_indel.vcf


#####################################################################################################

# ***IF EACH TOOLS' VARIANTS ARE TAKEN FOR TMB***

#sort TNsnv/${Patient}_TNsnv_completeFilter_SNP.txt Mutect2/${Patient}_Mutect2_completeFilter_SNP.txt TNscope/${Patient}_TNscope_completeFilter_SNP.txt | uniq > all_SNV.txt

#sort Mutect2/${Patient}_Mutect2_completeFilter_Ind.txt TNscope/${Patient}_TNscope_completeFilter_Ind.txt Strelka/${Patient}_Strelka_completeFilter_Ind.txt | uniq > all_Ind.txt

# Extract all SNPs and indels from the Mutect2 reference vcf and create a vcf ready for further annotation. See "8-coding_non-coding.sh" for further details.
#check=0
#for line in $(cat all_SNV.txt all_Ind.txt | sort -V)
#do
#if [ $check -eq 0 ]; then
#grep "#" Mutect2/${Patient}_Mutect2_preFilter.vcf > all_snv_indel.temp
#check=1
#fi
#grep -h "${line}" */${Patient}_*_preFilter.vcf >> all_snv_indel.temp
#done

#sed -e 's/Mutect2_//g' -e 's/Mutect2://g' -e 's/Strelka://g' -e 's/TNscope://g' -e 's/TNsnv://g' all_snv_indel.temp > all_snv_indel.vcf

#####################################################################################################

rm *.temp
rm annotation/*.temp
rm *.SNP
rm *.Ind








