#!/bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Script to run all other scripts included in the WGS analysis pipeline. This pipeline starts with the vcf files already acquired.

# Run from folder containing all patients folders. The patients folders must then contain the vcfs from tools used (Mutect2, Strelka, TNscope and TNsnv).
cd /mnt/raid5_hdd/ruelou01

# Some scripts already loop all patients, others loop over patients in this script.
# For more details, see each individual script.

sh 1.1-unfiltered_records.sh

for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//');
do

echo "\nExecuting patient ${Patient}"
Rscript 3-filtering_vcf_to_csv.R ${Patient}
sh 4-annotation.sh ${Patient}
sh 5-final_filtering.sh ${Patient}

done


sh 6-filtered_records.sh

for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//');
do

echo "\nExecuting patient ${Patient}"
sh 8-coding_non-coding.sh ${Patient}

done

sh 9-coding_non-coding_records.sh

sh 11-TMB.sh

# Scripts 2 and 7 (optional for Venn diagrams), 10, 12, 13 and 14 are ran in R (RStudio recommended).

# Comparing variant count before and after removal of dbSNP filter.
echo '"Patient_ID";"old_variant_count";"new_variant_count"'> old_dnSNP_vs_new_variant_count.csv

for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//');
do

old_variant_count=$(cat archive/Old_files_with_dbSNP_filtering_2022-07-06/NSLC_$Patient/common_2_snv_indel.vcf | wc -l)
new_variant_count=$(cat NSLC_$Patient/common_2_snv_indel.vcf | wc -l)

echo '"NSLC-'${Patient}'";"'${old_variant_count}'";"'${new_variant_count}'"' >> old_dnSNP_vs_new_variant_count.csv

done

# Converting vcf files into tsv files (for better readability in Excel).
for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//');
do

echo "\nExecuting patient ${Patient}"
grep -v '##' NSLC_${Patient}/NSLC_${Patient}_common_2_snv_indel.vcf | head -1 | sed -e 's/#//' -e 's/Mutect2_//g' > NSLC_${Patient}/NSLC_${Patient}_common_2_snv_indel.tsv
grep -v '#' NSLC_${Patient}/NSLC_${Patient}_common_2_snv_indel.vcf | sed 's/Mutect2_//g' >> NSLC_${Patient}/NSLC_${Patient}_common_2_snv_indel.tsv

done



