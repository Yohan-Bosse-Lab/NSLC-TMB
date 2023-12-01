#!/bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Calculating TMB of different sequencing regions and outputting results in a csv file.

cd /mnt/raid5_hdd/ruelou01

printf '"Patient_ID";"classic_WES_TMB";"complete_WES_TMB";"nc_WGS_TMB";"complete_WGS_TMB";"F1CDx_TMB"; "Illu500_TMB";"MSKI468_TMB";"NeoPlus_TMB";"OncoMineTMB_TMB";"QIAseq_TMB";"WES_mutation_count";"WGS_mutation_count"'> TMB.csv
printf "\n" >> TMB.csv

for Patient in $(ls vcf_Mutect2 | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//'); 
do

echo "\nCalculating Patient ${Patient} TMBs"

############
# WES TMBs #
############
echo "WES"
# Classic WES TMB (non-synonymous only)
# Exclusion of unmapped contigs variants (from GL000 contigs)
classic_WES=$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.exonic_variant_function | grep "nonsynonymous" | wc -l)
# the size of the exome is approximated to 32.8 Mb according to GATK (exons intervals available on GATK's public server https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
classic_WES_TMB=$(printf "%.3f\n" $(echo "${classic_WES} / 32.8" | bc -l))

# Complete WES TMB (all mutations except synonymous mutations)
# Exclusion of unmapped contigs variants (from GL000 contigs)
# Non-synonymous and framshifts count :
complete_WES_ns_fs=$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.exonic_variant_function | grep -w -v "synonymous" | wc -l)
# Splices sites count :
complete_WES_ss=$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.variant_function | grep "splicing" | wc -l)

WES_mutation_count=$(echo "${complete_WES_ns_fs}+${complete_WES_ss}" | bc -l)

# the size of the protein-coding region is approximated to 32.8 Mb according to GATK (CDS regions intervals available on GATK's public server https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
complete_WES_TMB=$(printf "%.3f\n" $(echo "${WES_mutation_count} / 32.8" | bc -l))

echo "WES mutation count: ${WES_mutation_count}"

############
# WGS TMBs #
############
echo "WGS"
# Non-coding WGS TMB (intronic and intergenic regions only)
# Exclusion of unmapped contigs variants (from GL000 contigs). Exonic mutations are excluded since we already have data from these mutations from the complete WES.
nc_WGS=$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.variant_function | grep -w -v "exonic" | wc -l)
# The size of the whole genome is approximated to 3000 Gb according to NCBI GRCh37 (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/#/st)
nc_WGS_TMB=$(printf "%.3f\n" $(echo "${nc_WGS} / 3000" | bc -l))

# Complete WGS TMB (all mutations excluding synonymous mutations)
# Exclusion of unmapped contigs variants (from GL000 contigs)
# The coding variant count is WES_mutation_count and the non-coding variant count is nc_WGS. At this stage, **synonymous mutations are excluded** since they are not taken into account in the complete WES.
WGS_mutation_count=$(printf "%.3f\n" $(echo "${WES_mutation_count}+${nc_WGS}" | bc -l))
# The size of the whole genome is approximated to 3 Gb according to NCBI GRCh37 (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/#/st)
complete_WGS_TMB=$(printf "%.3f\n" $(echo "${WGS_mutation_count} / 3000" | bc -l))

echo "WGS mutation count: ${WGS_mutation_count}"

###########
# NGS TMB #
###########

# Panel sizes were either provided by the manufacturer, online sources (namely https://www.researchgate.net/figure/Description-of-the-11-participating-diagnostic-NGS-panels_tbl1_340211075) or calculated with R 

# Genes panels NGS TMB 
for panel in $(ls NGS_panels_size/ | sed 's/_genes.txt//')
do

genes_count=0

for gene in $(cat NGS_panels_size/${panel}_genes.txt)
do
# Exclusion of unmapped contigs variants (from GL000 contigs) and search for each mutated gene count in available NGS panel genes. Includes all variants that are not synonymous in the annotation exonic file and also "splicing" variants (available in the non-exonic annotation file).
genes_count=$((genes_count+$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.exonic_variant_function | grep -w -v "synonymous" | grep -w "${gene}" | wc -l)+$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.variant_function | grep "splicing" | grep -w "${gene}" | wc -l)))
done

	if [ "$panel" = "F1CDx" ]
	then
	echo "${panel} gene count:${genes_count}"
	# the size of the F1CDx panel for TMB calculation is 0.8 Mb
	F1CDx_TMB=$(echo "${genes_count} / 0.8" | bc -l)
	fi

	if [ "$panel" = "Illu500" ]
	then
	echo "${panel} gene count:${genes_count}"
	# the size of the Illumina TruSeq Oncology 500 panel for TMB calculation is 1.33 Mb
	Illu500_TMB=$(echo "${genes_count} / 1.33" | bc -l)
	fi

	if [ "$panel" = "MSKI468" ]
	then
	echo "${panel} gene count:${genes_count}"
	# the size of the MSKImpact panel for TMB calculation is 1.14 Mb
	MSKI468_TMB=$(echo "${genes_count} / 1.14" | bc -l)
	fi

	if [ "$panel" = "NeoPlus" ]
	then
	echo "${panel} gene count:${genes_count}"
	# the size of the NeoPlus panel for TMB calculation is 1.1 Mb
	NeoPlus_TMB=$(echo "${genes_count} / 1.1" | bc -l)
	fi

	if [ "$panel" = "OncoMineTMB" ]
	then
	echo "${panel} gene count:${genes_count}"
	# the size of the OncoMineTMB panel for TMB calculation is 1.2 Mb
	OncoMineTMB_TMB=$(echo "${genes_count} / 1.2" | bc -l)
	fi

	if [ "$panel" = "QIAseq" ]
	then
	echo "${panel} gene count:${genes_count}"
	# the size of the QIAseq panel for TMB calculation is 1.33 Mb
	QIAseq_TMB=$(echo "${genes_count} / 1.33" | bc -l)
	fi

done

echo '"NSLC-'${Patient}'";"'${classic_WES_TMB}'";"'${complete_WES_TMB}'";"'${nc_WGS_TMB}'";"'${complete_WGS_TMB}'";"'${F1CDx_TMB}'";"'${Illu500_TMB}'";"'${MSKI468_TMB}'";"'${NeoPlus_TMB}'";"'${OncoMineTMB_TMB}'";"'${QIAseq_TMB}'";"'${WES_mutation_count}'";"'${WGS_mutation_count}'"' >> TMB.csv

done



