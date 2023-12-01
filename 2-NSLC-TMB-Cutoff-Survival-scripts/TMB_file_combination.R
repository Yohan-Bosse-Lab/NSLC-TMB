library(xlsx)
library(readxl)
library(dplyr)
library(openxlsx)

# Script for making a combined patients clinical and TMB data used for scripts 12 to 14.

# Importing several patient clinical characteristics.
#NSLC <- read_xlsx("../../NCI_NeverSmoker_n131_20181207_share.xlsx", range = "A1:V132") # Dating from 2018
NSLC <- read_xlsx("../NCI_NeverSmoker_n131_20210812_share_updated-124_125_fix.xlsx", range = "A1:W132") # Dating from 2021

# Importing patients ID (there are several identifications possible for each patient)
patients.id <- read_excel("../Table_ID-124_125_fix.xlsx", range = "A1:W125")[,c(1,8)]

# Importing updated tumor size, OS, pathological stage and comorbidities.
tumor_com_surv <- read_xlsx("TMB_tumor-size_comorbidities_LJR_2022-02-15.xlsx")
tumor_com_surv["comorbidities"] <- ifelse(tumor_com_surv$Asthma != "NA", "Asthma",
       ifelse(tumor_com_surv$Diabetes != "NA", "Diabetes",
              ifelse(tumor_com_surv$Emphysema != "NA", "Emphysema",
                     ifelse(tumor_com_surv$`HT - Hypertension` != "NA", "Hypertension",
                            ifelse(tumor_com_surv$`PH - Pulmonary hypertension` != "NA", "Pulmonary hypertension",
                                   ifelse(tumor_com_surv$COPD != "NA", "COPD",
                                          ifelse(tumor_com_surv$`Emphysema on CT scan` != "NA", "Emphysema", NA))))))) # Simplifying the comorbidities
tumor_com_surv <- tumor_com_surv[, -c(5:11)]

# Importing TMB values.
NSLC_TMB <- read.csv2("TMB.csv", stringsAsFactors = F)

# Updating clinical characteristics. already present in NSLC clinical file.
NSLC <- merge(NSLC, tumor_com_surv, by.x="ID", by.y = "NSujet")
NSLC$path.stage <- NSLC$Stage.final.vsa ; NSLC$Stage.final.vsa <- NULL # Updating pathological stages
colnames(NSLC)[24] <- "tumor.size (cm)"

# Merging clinical data with TMB data.
TMB_sample_ID <- merge(NSLC_TMB, patients.id, by.x="Patient_ID", by.y = "Subject")
NSLC_TMB <- merge(NSLC, TMB_sample_ID, by.x="normal_sample_ID", by.y = "Normal_Bam_ID")

# Reordering the variables.
NSLC_TMB <- NSLC_TMB[,c(27, 2, 1, 3:26, 28:length(colnames(NSLC_TMB)))]

# Changing TMB values to numeric (currently string)
NSLC_TMB[, 28:length(colnames(NSLC_TMB))] <- apply(NSLC_TMB[, 28:length(colnames(NSLC_TMB))], 2, function(x) as.numeric(x))
NSLC_TMB <- NSLC_TMB[order(NSLC_TMB$Patient_ID),]

# Writing new combined data table (containing patient clinical characteristics and TMB scores).
#write.xlsx(NSLC_TMB, "../NCI_NeverSmoker_n131_20210812_TMB.xlsx", row.names = FALSE)





### SAS data section

## ** OS and RpFS section.
# Explanation : the OS is already in the above dataset, but a more recent calculation was done and the unit was transferred into days (instead of months).

# Importing RpFS values
RpFS <- read_excel("4 - Survival analysis/Cohorte Never smookers_Calcul RFS et RpFS_2022-07-08_version LJ.xlsx", range = "A1:O94")
RpFS <- RpFS[,c("NSujet", "VitalStatus (0=vivant, 1=décédé)", "DDeces", "OS_calculé (jours)", "Indicateur RpFS (recidive=1, non=0)", "RpFS (jours)", "8TH Edition")]
RpFS <- RpFS[!RpFS$NSujet == 5804,] # Removing filtered patient (because of error in methodology)
colnames(RpFS) <- c("ID", "VitalStatus", "DDeces", "time_os", "RpFS_indicator", "time_RpFS", "pathological_stage")
RpFS$DDeces <- openxlsx::convertToDate(RpFS$DDeces)
RpFS$RpFS_indicator <- as.numeric(RpFS$RpFS_indicator)
RpFS$time_RpFS <- as.numeric(RpFS$time_RpFS)


NSLC_TMB <- NSLC_TMB[!NSLC_TMB$histology == 4,]
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$Patient_ID == "NSLC-0124",] # Removing filtered patient (because of error in methodology)
NSLC_TMB <- NSLC_TMB[!is.na(NSLC_TMB$normal_sample_ID),] # Removing NA row

NSLC_TMB <- merge(NSLC_TMB, RpFS, by.x="ID", by.y = "ID")

## Writing new combined data table for SAS survival analysis (containing patient clinical characteristics and TMB scores of the final 92 patients studied).
NSLC_TMB$'time.os (months)' <- NULL # The new OS data was integrated with the RpFS data.
NSLC_TMB$death <- NULL
NSLC_TMB$date.death <- NULL
NSLC_TMB$path.stage <- NULL
NSLC_TMB$date.death <- NULL
names(NSLC_TMB)[names(NSLC_TMB) == 'tumor.size (cm)'] <- 'tumor_size'
write.xlsx(NSLC_TMB, "4 - Survival analysis/SAS/NCI_NeverSmoker_n131_20210812_TMB_SAS.xlsx", rowNames = FALSE)
write.xlsx(NSLC_TMB[NSLC_TMB$histology == 3,], "../NCI_NeverSmoker_n131_20210812_TMB_adeno.xlsx", rowNames = FALSE)

