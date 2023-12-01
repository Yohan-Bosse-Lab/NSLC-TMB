#!/bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# This script compiles differents metrics (see first echo line) following the final filtering done in script 5.

cd /mnt/raid5_hdd/ruelou01

echo '"Patient_ID";"TNsnv_vcf";"Mutect2_vcf";"TNscope_vcf";"Strelka_vcf";"SNVs_common_to_3";"indels_common_to_3";"SNVs_common_to_2";"indels_common_to_2"'> Filtered_records.csv

for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//');
do

TNsnv_vcf=$(zcat NSLC_${Patient}/TNsnv/${Patient}_TNsnv_norm.vcf.gz | grep -v "#" | grep -v "GL000" | wc -l)
Mutect2_vcf=$(zcat NSLC_${Patient}/Mutect2/${Patient}_Mutect2_norm.vcf.gz | grep -v "#" | grep -v "GL000" | wc -l)
TNscope_vcf=$(zcat NSLC_${Patient}/TNscope/${Patient}_TNscope_norm.vcf.gz | grep -v "#" | grep -v "GL000" | wc -l)
Strelka_vcf=$(zcat NSLC_${Patient}/Strelka/${Patient}_Strelka_norm.indels.vcf.gz | grep -v "#" | grep -v "GL000" | wc -l)

# SNPs and indels counts of tools intersection for a each patient.
# Tools used for SNPs intersection : Mutect, Mutect2 and TNscope.
# Tools used for Indels intersection : Mutect2, TNscope and Strelka.

# SNPs section
grep -v "GL000" NSLC_${Patient}/TNsnv/${Patient}_TNsnv_completeFilter_SNP.txt | sort > NSLC_${Patient}/TNsnv_SNP_temp.txt
grep -v "GL000" NSLC_${Patient}/Mutect2/${Patient}_Mutect2_completeFilter_SNP.txt | sort > NSLC_${Patient}/Mutect2_SNP_temp.txt
grep -v "GL000" NSLC_${Patient}/TNscope/${Patient}_TNscope_completeFilter_SNP.txt | sort > NSLC_${Patient}/TNscope_SNP_temp.txt

# SNPs common to all 3 tools
snvs_common_3=$(comm -12 NSLC_${Patient}/TNsnv_SNP_temp.txt NSLC_${Patient}/Mutect2_SNP_temp.txt > NSLC_${Patient}/temp_SNP.txt; sort NSLC_${Patient}/temp_SNP.txt > NSLC_${Patient}/temp2_SNP.txt; comm -12 NSLC_${Patient}/temp2_SNP.txt NSLC_${Patient}/TNscope_SNP_temp.txt | wc -l)

# SNPs common to at least 2 tools
comm -12 NSLC_${Patient}/TNsnv_SNP_temp.txt NSLC_${Patient}/Mutect2_SNP_temp.txt | sort > NSLC_${Patient}/snvs_common_TNsnv_Mu2_temp.txt
comm -12 NSLC_${Patient}/TNsnv_SNP_temp.txt NSLC_${Patient}/TNscope_SNP_temp.txt | sort > NSLC_${Patient}/snvs_common_TNsnv_TNscope_temp.txt
comm -12 NSLC_${Patient}/Mutect2_SNP_temp.txt NSLC_${Patient}/TNscope_SNP_temp.txt | sort > NSLC_${Patient}/snvs_common_Mu2_TNscope_temp.txt

snvs_common_2=$(sort NSLC_${Patient}/snvs_common_TNsnv_Mu2_temp.txt NSLC_${Patient}/snvs_common_TNsnv_TNscope_temp.txt NSLC_${Patient}/snvs_common_Mu2_TNscope_temp.txt | uniq | wc -l)



# Indels section
grep -v "GL000" NSLC_${Patient}/Mutect2/${Patient}_Mutect2_completeFilter_Ind.txt | sort > NSLC_${Patient}/Mutect2_Ind_temp.txt
grep -v "GL000" NSLC_${Patient}/TNscope/${Patient}_TNscope_completeFilter_Ind.txt | sort > NSLC_${Patient}/TNscope_Ind_temp.txt
grep -v "GL000" NSLC_${Patient}/Strelka/${Patient}_Strelka_completeFilter_Ind.txt | sort > NSLC_${Patient}/Strelka_Ind_temp.txt

# Indels common to all 3 tools
indels_common_3=$(comm -12 NSLC_${Patient}/Mutect2_Ind_temp.txt NSLC_${Patient}/TNscope_Ind_temp.txt > NSLC_${Patient}/temp_Ind.txt; sort NSLC_${Patient}/temp_Ind.txt > NSLC_${Patient}/temp2_Ind.txt; comm -12 NSLC_${Patient}/temp2_Ind.txt NSLC_${Patient}/Strelka_Ind_temp.txt | wc -l)

# Indels common to at least 2 tools
comm -12 NSLC_${Patient}/Mutect2_Ind_temp.txt NSLC_${Patient}/TNscope_Ind_temp.txt | sort > NSLC_${Patient}/indels_common_Mu2_TNscope_temp.txt
comm -12 NSLC_${Patient}/Mutect2_Ind_temp.txt NSLC_${Patient}/Strelka_Ind_temp.txt | sort > NSLC_${Patient}/indels_common_Mu2_Strelka_temp.txt
comm -12 NSLC_${Patient}/TNscope_Ind_temp.txt NSLC_${Patient}/Strelka_Ind_temp.txt | sort > NSLC_${Patient}/indels_common_TNscope_Strelka_temp.txt

indels_common_2=$(sort NSLC_${Patient}/indels_common_Mu2_TNscope_temp.txt NSLC_${Patient}/indels_common_Mu2_Strelka_temp.txt NSLC_${Patient}/indels_common_TNscope_Strelka_temp.txt | uniq | wc -l)



echo "\nPatient $Patient SNPs and indels intersection count, reported in NSLC_${Patient}/Total_filtred_variants_${Patient}.txt.\n"
#echo "Total number of SNPs matching across Mutect, Mutect2 and TNscope for patient "${Patient}":" > NSLC_${Patient}/Total_filtered_variants_${Patient}.txt
#echo $snvs_common_3 >> NSLC_${Patient}/Total_filtered_variants_${Patient}.txt
#echo "\nTotal number of indels matching across Mutect, Mutect2 and TNscope for patient "${Patient}":" >> NSLC_${Patient}/Total_filtered_variants_${Patient}.txt
#echo $indels_common_3 >> NSLC_${Patient}/Total_filtered_variants_${Patient}.txt

echo "Total number of SNPs matching across at least 2 from Mutect, Mutect2 and/or TNscope for patient "${Patient}":" > NSLC_${Patient}/Total_filtered_variants_${Patient}.txt
echo $snvs_common_2 >> NSLC_${Patient}/Total_filtered_variants_${Patient}.txt
echo "\nTotal number of indels matching across at least 2 from Mutect, Mutect2 and TNscope for patient "${Patient}":" >> NSLC_${Patient}/Total_filtered_variants_${Patient}.txt
echo $indels_common_2 >> NSLC_${Patient}/Total_filtered_variants_${Patient}.txt

echo '"NSLC-'${Patient}'";"'${TNsnv_vcf}'";"'${Mutect2_vcf}'";"'${TNscope_vcf}'";"'${Strelka_vcf}'";"'${snvs_common_3}'";"'${indels_common_3}'";"'${snvs_common_2}'";"'${indels_common_2}'"' >> Filtered_records.csv

rm NSLC_${Patient}/*temp*.txt

done

