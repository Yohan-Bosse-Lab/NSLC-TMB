#!bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# DESC: This script prepares vcf files for further analyses. It compresses vcf files that aren't already compressed, does vairant normalization, rebuilds the indexes of normalized vcf files and annotates (adds easy identification of each variant) to simplify analyses.  
# This script is called in the script 1.1-unfiltered_records.sh. No need to use it by itself.
# For this script to work properly, change the paths of the vcf files corresponding to the different tools. Adjust the new directory paths if needed.

Patient=$1

if [ "$Patient" = "" ]; then
echo "Error: missing patient identifier argument."
exit
fi

# New directory paths
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}/SNPs
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}/Indels
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}/TNsnv
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}/Mutect2
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}/TNscope
mkdir /mnt/raid5_hdd/ruelou01/NSLC_${Patient}/Strelka

cd /mnt/raid5_hdd/ruelou01/NSLC_${Patient}

# Vcf files paths
Mutect="/mnt/raid5_hdd/ruelou01/vcf_TNsnv"
Mutect2="/mnt/raid5_hdd/ruelou01/vcf_Mutect2"
TNscope="/mnt/raid5_hdd/ruelou01/vcf_TNscope"
Strelka="/mnt/raid5_hdd/ruelou01/vcf_Strelka"

# File compression (.vcf to .vcf.gz), normalization, indexing section and annotation. Required for SNPs and indels stats.
# Annotation is done by changing the ID column (not used by the variant caller tools) to summarize each variant for simpler analyses.

# Gathering of normal and tumor ID to replace the them by corresponding "NORMAL" or "TUMOR" tags.
Tumor_ID=$(grep NSLC-${Patient} ../laval_wgs_124_paired_cram_124_125_fix-no_39.txt | awk '{print $4}')

Normal_ID=$(grep NSLC-${Patient} ../laval_wgs_124_paired_cram_124_125_fix-no_39.txt | awk '{print $5}')

echo "Patient $Patient TNsnv (Mutect) vcf preparation."
# Normalisation. Normal and tumor ID are changed to their corresponding tags and the archive is read beforehand.
sed "s/$Tumor_ID/TNsnv_TUMOR/" $Mutect/NSLC-${Patient}_TNsnv.vcf | sed "s/$Normal_ID/TNsnv_NORMAL/" | bgzip -c | bcftools norm -m -any -O z -o TNsnv/${Patient}_TNsnv_norm.vcf.gz
# Indexing
bcftools index -f TNsnv/${Patient}_TNsnv_norm.vcf.gz
# ID annotation 
bcftools annotate -O z -x 'ID' -I +'TNsnv:%CHROM:%POS:%REF:%ALT' TNsnv/${Patient}_TNsnv_norm.vcf.gz -o TNsnv/${Patient}_TNsnv_norm_ID.vcf.gz

echo "Patient $Patient TNhaplotyper2 (Mutect2) vcf preparation."
# Normalisation. Normal and tumor ID are changed to their corresponding tags and the archive is read beforehand.
sed "s/$Tumor_ID/Mutect2_TUMOR/" $Mutect2/NSLC-${Patient}_TNhaplotyper2.vcf | sed "s/$Normal_ID/Mutect2_NORMAL/" | bgzip -c | bcftools norm -m -any -O z -o Mutect2/${Patient}_Mutect2_norm.vcf.gz
# Indexing
bcftools index -f Mutect2/${Patient}_Mutect2_norm.vcf.gz
# ID annotation 
bcftools annotate -O z -x 'ID' -I +'Mutect2:%CHROM:%POS:%REF:%ALT' Mutect2/${Patient}_Mutect2_norm.vcf.gz -o Mutect2/${Patient}_Mutect2_norm_ID.vcf.gz

echo "Patient $Patient TNscope vcf preparation."
# Normalisation. Normal and tumor ID are changed to their corresponding tags and the archive is read beforehand.
sed "s/$Tumor_ID/TNscope_TUMOR/" $TNscope/NSLC-${Patient}_TNscope.vcf | sed "s/$Normal_ID/TNscope_NORMAL/" | bgzip -c | bcftools norm -m -any -O z -o TNscope/${Patient}_TNscope_norm.vcf.gz
# Indexing  
bcftools index -f TNscope/${Patient}_TNscope_norm.vcf.gz
# ID annotation 
bcftools annotate -O z -x 'ID' -I +'TNscope:%CHROM:%POS:%REF:%ALT' TNscope/${Patient}_TNscope_norm.vcf.gz -o TNscope/${Patient}_TNscope_norm_ID.vcf.gz

echo "Patient $Patient Strelka vcf preparation."
# Normalisation. The archive is read beforehand. "NORMAL" and "TUMOR" tags are already present in Strelka vcf files.
zcat $Strelka/NSLC-${Patient}_strelka/results/variants/somatic.indels.vcf.gz | sed "s/TUMOR/Strelka_TUMOR/" | sed "s/NORMAL/Strelka_NORMAL/" | bgzip -c | bcftools norm -m -any -O z -o Strelka/${Patient}_Strelka_norm.indels.vcf.gz
# Indexing
bcftools index -f Strelka/${Patient}_Strelka_norm.indels.vcf.gz
# ID annotation 
bcftools annotate -O z -x 'ID' -I +'Strelka:%CHROM:%POS:%REF:%ALT' Strelka/${Patient}_Strelka_norm.indels.vcf.gz -o Strelka/${Patient}_Strelka_norm_ID.indels.vcf.gz
echo "Preparation done\n"

# SNPs and indels records for Venn diagram use
echo "Processing SNPs and indels records, reported in directories NSLC-$Patient/SNPs/ and NSLC-$Patient/Indels/.\n"

# SNPs count for TNsnv
bcftools view TNsnv/${Patient}_TNsnv_norm_ID.vcf.gz | cut -f3 | grep -v "^##" | grep -v "^ID" | awk 'BEGIN {FS = ":"} ; {if (length($4)==1 && length($5)==1){print > "SNPs/TNsnv_SNP.txt"}}'

# SNPs and indels count for TNhaplotyper2 (Mutect2)
bcftools view Mutect2/${Patient}_Mutect2_norm_ID.vcf.gz | cut -f3 | grep -v "^##" | grep -v "^ID" | awk 'BEGIN {FS = ":"} ; {if (length($4)==1 && length($5)==1){print > "SNPs/Mutect2_SNP.txt"}; if (length($4)!=length($5)){print > "Indels/Mutect2_Ind.txt"}}'

# SNPs and indels count for TNscope
bcftools view TNscope/${Patient}_TNscope_norm_ID.vcf.gz | cut -f3 | grep -v "^##" | grep -v "^ID" | awk 'BEGIN {FS = ":"} ; {if (length($4)==1 && length($5)==1){print > "SNPs/TNscope_SNP.txt"}; if (length($4)!=length($5)){print > "Indels/TNscope_Ind.txt"}}'

# Indels count for Strelka
bcftools view Strelka/${Patient}_Strelka_norm_ID.indels.vcf.gz | cut -f3 | grep -v "^##" | grep -v "^ID" | awk 'BEGIN {FS = ":"} ; {if (length($4)!=length($5)){print > "Indels/Strelka_Ind.txt"}}'



