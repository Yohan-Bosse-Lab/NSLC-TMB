#!/bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# DESC: This script reads all the variants called by TNsnv, TNhaplotpyer2 (Mutect2) TNscope and Strelka. They are then documented in the file "Unfiltered_records.csv". 
# This script also calls the complementary script 1.2-vcf_prep.sh. No need to use the latter by itself.

cd /mnt/raid5_hdd/ruelou01

# Preparing the file to document all unfiltered variant records. The filtering takes place in scripts 3, 4 and 5. The filters include read depth, variant read count and variant allele frequency thresholds.
echo '"Patient_ID";"TNsnv_vcf";"Mutect2_vcf";"TNscope_vcf";"Strelka_vcf";"SNVs_common";"indels_common"'> Unfiltered_records.csv

# Going through the entire patients list. I take Mutect2 file names as a reference, but could be another reference list.
for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//'); do

# Run the vcf preparation script to have streamlined vcfs files for the current patient.
sh 1.2-vcf_prep.sh $Patient

# Import streamlined vcf file for the current patient
TNsnv_vcf=$(zcat NSLC_${Patient}/TNsnv/${Patient}_TNsnv_norm.vcf.gz | grep -v "#" | wc -l)
Mutect2_vcf=$(zcat NSLC_${Patient}/Mutect2/${Patient}_Mutect2_norm.vcf.gz | grep -v "#" | wc -l)
TNscope_vcf=$(zcat NSLC_${Patient}/TNscope/${Patient}_TNscope_norm.vcf.gz | grep -v "#" | wc -l)
Strelka_vcf=$(zcat NSLC_${Patient}/Strelka/${Patient}_Strelka_norm.indels.vcf.gz | grep -v "#" | wc -l)

# SNPs and indels counts of tools intersection for a specific patient.
# Tools used for SNPs intersection : Mutect, Mutect2 and TNscope.
# Tools used for Indels intersection : Mutect2, TNscope and Strelka.
snvs_common=$(bcftools isec -n-3 NSLC_${Patient}/TNsnv/${Patient}_TNsnv_norm.vcf.gz NSLC_${Patient}/Mutect2/${Patient}_Mutect2_norm.vcf.gz NSLC_${Patient}/TNscope/${Patient}_TNscope_norm.vcf.gz | awk '{if ($5=="111"){sum+=1}}; END{print sum}')

indels_common=$(bcftools isec -n-3 NSLC_${Patient}/Mutect2/${Patient}_Mutect2_norm.vcf.gz NSLC_${Patient}/TNscope/${Patient}_TNscope_norm.vcf.gz NSLC_${Patient}/Strelka/${Patient}_Strelka_norm.indels.vcf.gz | awk '{if ($5=="111"){sum+=1}}; END{print sum}')

# Reporting common SNPs and common indels across variant calling tools.
echo "\nPatient $Patient SNPs and indels intersection count, reported in NSLC_${Patient}/Total_unfiltred_variants_${Patient}.txt.\n"

echo "Total number of SNPs matching across Mutect, Mutect2 and TNscope for patient "${Patient}":" > NSLC_${Patient}/Total_unfiltred_variants_${Patient}.txt
echo $snvs_common >> NSLC_${Patient}/Total_unfiltred_variants_${Patient}.txt

echo "\nTotal number of indels matching across Mutect, Mutect2 and TNscope for patient "${Patient}":" >> NSLC_${Patient}/Total_unfiltred_variants_${Patient}.txt
echo $indels_common >> NSLC_${Patient}/Total_unfiltred_variants_${Patient}.txt

# Final reporting of variants across all tools, common SNPs and common indels.
echo '"NSLC-'${Patient}'";"'${TNsnv_vcf}'";"'${Mutect2_vcf}'";"'${TNscope_vcf}'";"'${Strelka_vcf}'";"'${snvs_common}'";"'${indels_common}'"' >> Unfiltered_records.csv

done

