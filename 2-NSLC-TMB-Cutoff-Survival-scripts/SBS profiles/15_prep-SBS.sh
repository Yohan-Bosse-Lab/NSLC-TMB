#!bin/bash

# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Preparing files to run SBS profiles annotation on all patients.

cd /mnt/raid5_hdd/ruelou01/
mkdir 92_patients_NSLC_filtered_VCFS

for Patient in $(cat 92_patients_list.csv | tail -n +2 | sed 's/"//g');
do

echo $Patient
cp NSLC_${Patient}/NSLC_${Patient}_common_2_snv_indel.vcf 92_patients_NSLC_filtered_VCFS/NSLC_${Patient}.vcf

done
