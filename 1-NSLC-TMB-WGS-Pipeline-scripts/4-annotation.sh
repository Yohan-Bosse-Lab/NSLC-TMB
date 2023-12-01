#!/bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Annotation of vcf files by filtering them against dbSNP151, 1000 genome phase 3 v5, ExAC v0.3.1 and gnomAD v2.1.1 databases

# Annovar database names
#dbSNP151 = snp151
#1000 genome = 1000g2015aug_all
#ExAC v0.3.1 = exac03
#gnomAD = gnomad211_genome

Patient=$1

if [ "$Patient" = "" ]; then
echo "Error: missing patient identifier argument."
exit
fi

cd /mnt/raid5_hdd/ruelou01/NSLC_${Patient}

mkdir annotation

# TNsnv, Mutect2 and Tnscope vcf files compression, required for merging all files
for tool in TNsnv Mutect2 TNscope
do

bgzip -c ${tool}/${Patient}_${tool}_preFilter.vcf > ${tool}/${Patient}_${tool}_preFilter.vcf.gz
tabix ${tool}/${Patient}_${tool}_preFilter.vcf.gz

if [ "$tool" = "TNsnv" ]
then
# vcfout.list is the list of files to merge
ls -R ${tool}/*preFilter.vcf.gz > vcfout.list
else
ls -R ${tool}/*preFilter.vcf.gz >> vcfout.list
fi

done

# Strelka file modification and compression
# Strelka vcf files are missing the GT (genotype) FORMAT entry in order for Annovar to work.
# The first few steps are to add the GT entry in format header and in samples.
# Default genotypes 0/0 and 0/1 are applied for normal and tumor samples respectively. These values do not affect Annovar's output. 
strelka_formal="Strelka/${Patient}_Strelka_preFilter.vcf"
strelka_formal_mod="Strelka/${Patient}_Strelka_preFilter_GT.vcf"
first_format_num=$(grep -n -m 1 '##FORMAT' "$strelka_formal" | cut -d : -f 1)
sed "$first_format_num"'i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' "$strelka_formal" > "$strelka_formal_mod"
sed -ri 's|(DP:)|GT:\1|g' "$strelka_formal_mod"
sed -ri 's|(:BCN50\t)|\10/0:|g' "$strelka_formal_mod"
sed -ri 's|(:BCN50\t[^\t]*\t)|\10/1:|g' "$strelka_formal_mod"

# Remaking the Strelka vcf archive.
bgzip -c Strelka/${Patient}_Strelka_preFilter_GT.vcf > Strelka/${Patient}_Strelka_preFilter_GT.vcf.gz
tabix Strelka/${Patient}_Strelka_preFilter_GT.vcf.gz

echo ${strelka_formal_mod}.gz >> vcfout.list



# Merge of compressed files from TNsnv, Mutect2, TNscope and Strelka
bcftools merge -m id --file-list vcfout.list -O v -o annotation/preFilterMerge.vcf

# Making of inverted REF and ALT annovar input file from merged file. This is done if the "ALT" allele is the one found in databases, instead of the "REF" one.
grep '#' annotation/preFilterMerge.vcf > annotation/preFilterMerge_inverted.vcf
grep -v '#' annotation/preFilterMerge.vcf | awk 'BEGIN{OFS="\t"} {t = $4; $4 = $5; $5 = t; print}' >> annotation/preFilterMerge_inverted.vcf



# Annotation of formal (non-inverted) and inverted merged vcf file
/mnt/raid5_hdd/annovar/table_annovar.pl annotation/preFilterMerge.vcf /mnt/raid5_hdd/annovar/humandb/ -buildver hg19 -remove -protocol snp151,1000g2015aug_all,exac03,gnomad211_genome -operation f,f,f,f -vcfinput -out annotation/formal_annotation & 
/mnt/raid5_hdd/annovar/table_annovar.pl annotation/preFilterMerge_inverted.vcf /mnt/raid5_hdd/annovar/humandb/ -buildver hg19 -remove -protocol snp151,1000g2015aug_all,exac03,gnomad211_genome -operation f,f,f,f -vcfinput -out annotation/inverted_annotation & 
wait

rm vcfout.list


