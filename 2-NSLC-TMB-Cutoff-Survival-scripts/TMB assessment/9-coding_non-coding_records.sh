#!bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Merging of all coding and non-coding vcf files for all patients. This will report all significant variants across all patients, wether they are coding or non-coding.
# 2 directories are outputted : 
# - coding_non-coding_fullRecords/ for exact common variants between two patients -> comparison by using full variants records, including gene name, position, description, ...
# - coding_non-coding_genesRecords/ for common variants between two patients -> comparison by using chromosome and gene recrds. This informs of common mutated genes between two patients without any variant specificity 

cd /mnt/raid5_hdd/ruelou01

mkdir coding_non-coding_fullRecords
mkdir coding_non-coding_genesRecords

# Preparing sorted full records and genes records for all patients. Required for the next loop where a patient and the following patient are analysed. The following patient file must already be sorted, hence this loop.
# Unmapped contigs variants (GL000 contigs) are excluded from patients comparisons.
for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//'); 
do

sed -e 1d NSLC_${Patient}/annotation/coding_non-coding_${Patient}.csv | sort | grep -v "GL000" > NSLC_${Patient}/annotation/coding_non-coding_${Patient}.temp

awk 'BEGIN{FS=";"} {print $1, $7}' NSLC_${Patient}/annotation/coding_non-coding_${Patient}.temp | sort > NSLC_${Patient}/annotation/coding_non-coding_${Patient}_genes.temp

done


# Table creation of coding and non-coding variant counts for all patients
printf '"Patient_ID";"filtered_SNVs_common";"filtered_indels_common";"coding_mutations";"noncoding_mutations"' > coding_non-coding_records.csv
printf "\n" >> coding_non-coding_records.csv

# List of patients in order to get the next patient following the current one. Required to find overlapping variants.
ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//' > patient.list 

for Patient in $(ls vcf_Mutect2/ | grep -v ".vcf.idx" | awk 'BEGIN{FS="_"}{print $1}' | sed 's/NSLC-//')
do

for NextPatient in $(sed -n "/${Patient}/,\$p" patient.list | tail -n+2)
do

echo Comparing $Patient with $NextPatient.
comm -12 NSLC_${Patient}/annotation/coding_non-coding_${Patient}.temp NSLC_${NextPatient}/annotation/coding_non-coding_${NextPatient}.temp | sort -V > coding_non-coding_fullRecords/${Patient}_${NextPatient}

comm -12 NSLC_${Patient}/annotation/coding_non-coding_${Patient}_genes.temp NSLC_${NextPatient}/annotation/coding_non-coding_${NextPatient}_genes.temp | sort -V > coding_non-coding_genesRecords/${Patient}_${NextPatient}

done

# Count of coding and non-coding variants.
# For the coding variants, we do not count mutations from unmapped contigs or synonymous mutations and count the rest.
coding_count=$(($(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.exonic_variant_function | grep -w -v "synonymous" | grep -w "${gene}" | wc -l)+$(grep -v "GL000" NSLC_${Patient}/annotation/coding_non-coding_${Patient}.variant_function | grep "splicing" | wc -l)))

# For the non-coding variants, we do not count any coding variant present in the ".variant_function" file. They are annoted as "exonic" or "splicing" in the fist column.
noncoding_count=$(($(grep -v "splicing" < NSLC_${Patient}/annotation/coding_non-coding_${Patient}.variant_function | wc -l)-$(wc -l < NSLC_${Patient}/annotation/coding_non-coding_${Patient}.exonic_variant_function)))

# Reporting counts while evaluating variants common to all 3 tools.
#grep "NSLC-${Patient}" Filtered_records.csv | sed 's/\"//g' | awk -v cc="${coding_count}" -v ncc="${noncoding_count}" 'BEGIN{FS=";";OFS="\";\""} {print "\""$1,$6,$7,cc,ncc"\""}' >> coding_non-coding_records.csv

# Reporting counts while evaluating variants common to at least 2 tools.
grep "NSLC-${Patient}" Filtered_records.csv | sed 's/\"//g' | awk -v cc="${coding_count}" -v ncc="${noncoding_count}" 'BEGIN{FS=";";OFS="\";\""} {print "\""$1,$8,$9,cc,ncc"\""}' >> coding_non-coding_records.csv

done

rm NSLC_*/annotation/coding_non-coding_*.temp
rm patient.list

# Removing empty comparisons.
find coding_non-coding_fullRecords -size  0 -delete
find coding_non-coding_genesRecords -size  0 -delete
echo "\nNull comparisons (size of 0) have been deleted. The remaining comparisons are present in coding_non-coding_fullRecords/ and coding_non-coding_genesRecords/ directories."



