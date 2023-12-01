
library(xlsx)
library(data.table)
library(tidyverse)
library(glue)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
########################### Clinicopathological data ###########################
################################################################################

# Updated PFS and OS data (2023-04-20)
survival_data <- as.data.table(read.xlsx("Data/Never smokers_pfs_os_mise à jour 2022_2023-04-20_vsa.xlsx", sheetIndex = 2, endRow = 93, colIndex = c(1:5)))

# Clinicopathological data
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))

# Updating the clinicopathological data with new PFS and OS data
clin_data[, c("VitalStatus", "DDeces", "time_os", "time_RpFS", "RpFS_indicator") := NULL] # Removing old variables
clin_data <- merge(clin_data, survival_data, by.x = "Patient_ID", by.y = "Patient_ID") # Adding new variables
setcolorder(clin_data, names(clin_data)[c(1:14, 46:49, 15:45)]) # Reordering columns

# Adjusting PFS: if patient is dead (VitalStatus == 1), then PFS indicator should be 1.
clin_data[VitalStatus == 1 & PFS_indicator == 0, PFS_indicator := 1]

# Keeping only adenocarcinomas
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# TMB data for whole genome and specific regions
TMB <- as.data.table(read.xlsx("Data/VEP_82_NSLC_TMB.xlsx", sheetIndex = 1)) # created in the script 1-82_NSLC_vcf_TMB.r

# Merged dataset for survival models
surv.dt <- merge.data.table(clin_data, TMB, by = "Patient_ID")
surv.dt <- surv.dt %>% mutate(TMB_high_low = case_when(complete_WGS_TMB >= 1.70 ~ "High",
                                                       complete_WGS_TMB < 1.70 ~ "Low"),
                              pathological_stage_refactor = case_when(pathological_stage %in% c("1A1","1A2","1A3","1B") ~ "I",
                                                                      pathological_stage %in% c("2A", "2B") ~ "II",
                                                                      pathological_stage %in% c("3A", "3B", "4A") ~ "III & IV"),
                              pfs_two = case_when(PFS_time >= 730 ~ ">=2",
                                                  PFS_time < 730 ~ "<2"),
                              pfs_five = case_when(PFS_time >= 1825 ~ ">=5",
                                                   PFS_time < 1825 ~ "<5"),
                              os_two = case_when(OS_time >= 730 ~ ">=2",
                                                 OS_time < 730 ~ "<2"),
                              os_five = case_when(OS_time >= 1825 ~ ">=5",
                                                  OS_time < 1825 ~ "<5"))

surv.dt <- surv.dt %>% relocate(pathological_stage_refactor, .after = pathological_stage)

# Changing comorbidities NA to None
surv.dt[, comorbidities:=replace_na(comorbidities, "None")]

# Replacing dots in variable name for underscore for better portability
colnames(surv.dt) <- gsub("[.]", "_", colnames(surv.dt))

# Exporting the final clinicopathological dataset
write.xlsx(surv.dt, "Data/n82_NSLC_surv_data.xlsx", row.names = FALSE)

################################################################################
################################ Variants data #################################
################################################################################

# Fetching and merging all patients variants
VEP_data_all_patients <- rbindlist(lapply(Patients_list, function(x) fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{x}.csv"))), use.names = TRUE)
# Adding patient ID (for easier access than the excisting "ID" field)
VEP_data_all_patients[, Patient_ID := as.factor(gsub(":.*", "", ID))]
VEP_data_all_patients <- VEP_data_all_patients %>% relocate(Patient_ID, .after = ID)

# Exporting the final variants dataset
fwrite(VEP_data_all_patients, "Data/n82_NSLC_variants_data.csv", row.names = FALSE)
