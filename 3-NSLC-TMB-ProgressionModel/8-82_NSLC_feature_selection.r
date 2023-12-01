
# Question: Faut-il aussi inclure les patients qui n'ont pas eu d'evenemement avant 2 ans, afin d'aider le modele ? faire la difference chez les patients?
# Reponse: Non, car on ne peut pas tirer de conclusion sur ces patients (s'ils ont eu une evenement ou non ? l'interieur de 2 ans), car le suivi est trop court (< 2 ans).

library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(survival)
library(pls)
library(progress)
library(caret)
library(mlbench)
library(MXM)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################## Data import and formatting ############################
################################################################################

# Window size used to compute sliding windows
window_size <- 5000
Mb_length <-  window_size/1e6

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/n82_NSLC_surv_data.xlsx", sheetIndex = 1))

# Event of interest: progression event (relapse or death) before 2 years post-surgery
clin_data[, PFS_indicator_adjusted := ifelse(PFS_indicator == 1 & pfs_two == "<2", 1, 0)]
# Keeping only patients with an event before 2 years or at least 2 years of follow-up
studied_patients_data <- clin_data[!(PFS_indicator == 0 & pfs_two == "<2"), ]
patient_list <- gsub("NSLC-", "", studied_patients_data$Patient_ID)
studied_patients_data[, Patient_ID := patient_list]

## Merging all chromosome sliding windows densities together
chrom_name <- c(1:22, 'X')

# Reading density data of chromosomes with a specific window size
chrom_frequency <- rbindlist(lapply(chrom_name,
                                    function(x) fread(glue("Data/Sliding windows and frequency - individual patients/NSLC_{x}_sw_density_{Mb_length}Mb.csv"))),
                             use.names = TRUE)

# Keep only sliding windows with studied patients
chrom_frequency <- chrom_frequency %>%
  filter(str_detect(Patients_ID, paste(patient_list, collapse = "|")))

# Subtracting the variant frequency of sliding windows containing excluded patients
excluded_patients <- clin_data[!(gsub("NSLC-", "", Patient_ID) %in% patient_list), gsub("NSLC-", "", Patient_ID)]
excluded_patients.index <- str_count(chrom_frequency$Patients_ID, paste(excluded_patients, collapse = "|"))
chrom_frequency[, variant_frequency_excluded := variant_frequency - excluded_patients.index]
chrom_frequency[, variant_frequency := NULL]

################################################################################
############################## Preparing data ##################################
################################################################################

# # Separating sliding windows by each patient
# split_patient_data <- chrom_frequency %>% separate_longer_delim(Patients_ID, delim=";")
# setDT(split_patient_data)
# rm(chrom_frequency)
# 
# # Adding a window ID for each entry to replaces other fields. The syntax is chromosome.window_start.window_end.
# split_patient_data[, window_ID := paste(CHROM, window_start, window_end, sep="."), by = list(CHROM, window_start)]
# split_patient_data_bpk <- split_patient_data
# # Removing repetitive columns to free some space
# split_patient_data[, `:=`(window_start = NULL, window_end = NULL, CHROM = NULL, variant_frequency = NULL, variant_frequency_excluded = NULL)]
# 
# # Counting variants from each patient for each window
# split_patient_data[, frequency := .N, by=list(window_ID, Patients_ID)]
# 
# # Reordering columns
# split_patient_data <- split_patient_data %>% relocate(window_ID, .before = Patients_ID)
# 
# # Removing duplicated rows after counting frequencies
# split_patient_data.duplicates <- !duplicated(split_patient_data)
# split_patient_data <- split_patient_data[split_patient_data.duplicates]
# 
# # Removing entries with a variant frequency of 0. They are identified by an empty Patients_ID field.
# split_patient_data <- split_patient_data[Patients_ID != ""]
# 
# # Changing patient ID column name for better concordance between future datasets
# colnames(split_patient_data)[2] <- "Patient_ID"
# 
# # Checkpoint to free memory
# fwrite(split_patient_data, glue("Data/Sliding windows and frequency - individual patients/split_variant_data_{Mb_length}Mb.csv"))

# # Window size used to compute sliding windows
# window_size <- 5000
# Mb_length <-  window_size/1e6
#
# # Opening checkpoint file
# split_patient_data <- fread(glue("Data/Sliding windows and frequency - individual patients/split_variant_data_{Mb_length}Mb.csv"), colClasses = c('Patient_ID'='character'))
#
# # Generating the table used for partial least square (PLS) regression
# PLS_data <- as.data.table(split_patient_data %>%
#   group_by(window_ID) %>%
#   pivot_wider(
#     names_from = window_ID,
#     values_from = frequency))
#
# # Freeing memory
# rm(split_patient_data)
#
# PLS_data <- PLS_data[Patient_ID %in% gsub("NSLC-", "", patient_list),]
# PLS_data <- PLS_data[order(Patient_ID)]
#
# # Checkpoint to free memory
# fwrite(PLS_data, glue("Data/Sliding windows and frequency - individual patients/PLS_data_{Mb_length}Mb.csv"))

# # Window size used to compute sliding windows
# window_size <- 5000
# Mb_length <-  window_size/1e6
# # Opening checkpoint file
# PLS_data <- fread(glue("Data/Sliding windows and frequency - individual patients/PLS_data_{Mb_length}Mb.csv"), colClasses = c('Patient_ID'='character'))
#
# # Fast function to replace NAs in a data.table
# replace_na_dt = function(DT) {
#   for (j in seq_along(DT))
#     set(DT,which(is.na(DT[[j]])),j,0)
# }
#
# replace_na_dt(PLS_data)
#
# # Adding the PFS indicator to the PLS data
# PLS_data <- merge(studied_patients_data[, .(Patient_ID, PFS_indicator_adjusted, PFS_time)], PLS_data, by="Patient_ID")
#
# # Final data for PLS use
# fwrite(PLS_data, glue("Data/Sliding windows and frequency - individual patients/PLS_data_{Mb_length}Mb.csv"))

################################################################################
############################# ALL EVENT PATIENTS ###############################
################################################################################

# ## Preparing data
# data.read <- fread(glue("Data/Sliding windows and frequency - individual patients/PLS_data_{Mb_length}Mb.csv"), colClasses = c('Patient_ID'='character'))
# 
# ## Prefiltering the windows to reduce the dimensionality of the regression
# # Keeping 50k windows with the highest frequency, without filtering for the PFS indicator
# PLS_data_filter <- studied_patients_data[, .(Patient_ID, PFS_indicator_adjusted)]
# PLS_data_filter[, Patient_ID := gsub("NSLC-", "", Patient_ID)]
# most_mutated_sw.id <- chrom_frequency[, order(-variant_frequency_excluded)][1:50000]
# most_mutated_sw.data <- chrom_frequency[most_mutated_sw.id]
# 
# pb <- progress_bar$new(format = "[:bar] :percent",
#                        total = nrow(most_mutated_sw.data),
#                        complete = "=",   # Completion bar character
#                        incomplete = "-", # Incomplete bar character
#                        current = ">",    # Current bar character
#                        clear = FALSE,    # If TRUE, clears the bar when finish
#                        width = 100)      # Width of the progress bar
# 
# for(window in 1:nrow(most_mutated_sw.data)) {
#   pb$tick()
#   
#   if(most_mutated_sw.data[window]$variant_frequency == 0) next
#   
#   window_name <- paste(c(most_mutated_sw.data[window]$CHROM, most_mutated_sw.data[window]$window_start, most_mutated_sw.data[window]$window_end), collapse=".")
#   window_patients <- as.data.table(plyr::count(unlist(str_split(most_mutated_sw.data[window, .(Patients_ID)], ";"))))
#   colnames(window_patients) <- c("Patient_ID", window_name)
#   if(window %% 1000 == 0) {
#     start.time <- Sys.time()
#     PLS_data_filter[, glue(window_name) := merge(x=PLS_data_filter[, .(Patient_ID)], y=window_patients, by = "Patient_ID", all.x = TRUE)[[glue(window_name)]]]
#     end.time <- Sys.time()
#     time.taken <- end.time - start.time
#     print(time.taken)
#   } else {
#     PLS_data_filter[, glue(window_name) := merge(x=PLS_data_filter[, .(Patient_ID)], y=window_patients, by = "Patient_ID", all.x = TRUE)[[glue(window_name)]]]
#   }
# }
# 
# colnames(PLS_data_filter) <- as.character(gsub("\\.", ":", colnames(PLS_data_filter)))
# 
# # Replacing NAs mutation frequency with 0
# replace_na_dt = function(DT) {
#   for (j in seq_along(DT))
#     set(DT,which(is.na(DT[[j]])),j,0)
# }
# 
# replace_na_dt(PLS_data_filter)
# 
# fwrite(PLS_data_filter, glue("Data/Sliding windows and frequency - individual patients/50k_most_mutated_sw_{Mb_length}Mb.csv"))


################################################################################
# Import prefiltered data: the 50k highest mutation frequency sliding windows
data.read.filter <- fread(glue("Data/Sliding windows and frequency - individual patients/50k_most_mutated_sw_{Mb_length}Mb.csv"), 
                          colClasses = c('Patient_ID'='character'), stringsAsFactors = FALSE)


################################################################################
## Using caret for feature selection
# Data is not split between train and test because of the low positive sample size (n=10, control n=63). This is a major limitation.

x <- data.read.filter[, 3:ncol(data.read.filter)]

normalization <- preProcess(x)
x <- predict(normalization, x)
subsets <- seq(1000, 50000, 1000)

y <- as.factor(data.read.filter$PFS_indicator_adjusted)



## Using cross-validation
# define the control using a random forest selection function
control.cv <- rfeControl(functions=rfFuncs, method="cv", number = 10)
# run the RFE algorithm
results.cv <- rfe(x, y, sizes = subsets, rfeControl=control.cv)
# summarize the results
print(results.cv)
# list the chosen features
results.cv.final.sw <- str_sort(predictors(results.cv), numeric = TRUE)
fwrite(list(results.cv.final.sw), file=glue("Results/Feature Selection/cv_{Mb_length}Mb_final_features_all_patients.csv"))

# plot the results
ggplot(results.cv) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  theme_classic()

ggsave(file=glue("Results/Feature Selection/cv_{Mb_length}Mb_50k_features_all_patients.png"), width=6, height=6, dpi=300)



## Using boostrap
control <- rfeControl(functions = rfFuncs, method = "boot", verbose = TRUE,
                      returnResamp = "final", number = 10)

if ( require("parallel", quietly = TRUE, warn.conflicts = FALSE) ) {
  control$workers <- parallel:::detectCores()
  control$computeFunction <- mclapply
  control$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
}

results.boot <- rfe(x, y, sizes = subsets, rfeControl=control)
# summarize the results
print(results.cv)
# list the chosen features
print(results.boot)
results.boot.final.sw <- str_sort(predictors(results.boot), numeric = TRUE)
fwrite(list(results.boot.final.sw), file=glue("Results/Feature Selection/boot_{Mb_length}Mb_final_features_all_patients.csv"))

# plot the results
ggplot(results.boot) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  theme_classic()

ggsave(file=glue("Results/Feature Selection/boot_{Mb_length}Mb_50k_features_all_patients.png"), width=6, height=6, dpi=300)



################################################################################
## Using MXM for feature selection : https://cran.r-project.org/web/packages/MXM/vignettes/MMPC_tutorial.html
target <- as.factor(data.read.filter$PFS_indicator_adjusted)

dataset <- data.read.filter[, 3:ncol(data.read.filter)]
normalization <- preProcess(dataset)
dataset <- predict(normalization, dataset)

mmpc <- MMPC(target = target,
             dataset = dataset,
             max_k = 1000,
             threshold = 0.05,
             test = "testIndLogistic",
             ini = NULL,
             hash = TRUE,
             hashObject = NULL,
             ncore = 4,
             backward = TRUE)

summary(mmpc)
mmpc@selectedVarsOrder

dataset[,39817]

supervised.pca(target, dataset, indices = 1:5000)




################################################################################
################################## ARCHIVES ####################################
################################################################################

# Choice of supervised PCA : https://stats.stackexchange.com/questions/108614/regression-in-pn-setting-how-to-choose-regularization-method-lasso-pls-pc

# Using Superpc for supervised PCA feature selection
# Preparing the data for Supervised PCA

options(scipen=999)

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/n82_NSLC_surv_data.xlsx", sheetIndex = 1))

# Event of interest: progression event (relapse or death) before 2 years post-surgery
clin_data[, PFS_indicator_adjusted := ifelse(PFS_indicator == 1 & pfs_two == "<2", 1, 0)]
# Keeping only patients with an event before 2 years or at least 2 years of follow-up
studied_patients_data <- clin_data[!(PFS_indicator == 0 & pfs_two == "<2"), ]
patient_list <- gsub("NSLC-", "", studied_patients_data$Patient_ID)
studied_patients_data[, Patient_ID := patient_list]

# Adding PFS_time to the filtered data
data.read.filter <- merge(data.read.filter, studied_patients_data[, .(Patient_ID, PFS_time, TMB_high_low)], by="Patient_ID")
data.read.filter <- data.read.filter %>% relocate(c(PFS_time, TMB_high_low), .after = PFS_indicator_adjusted) # Relocating PFS_time after PFS_indicator

# Preparing the data for Supervised PCA
x <- data.read.filter[, 5:ncol(data.read.filter)]
normalization <- preProcess(x)
x <- predict(normalization, x)
x <- t(as.matrix(x))

y <- rnorm(data.read.filter$PFS_time)
censoring.status <- data.read.filter$PFS_indicator_adjusted
featurenames <- colnames(data.read.filter)[5:ncol(data.read.filter)]

data <- list(x = x,
             y = y,
             censoring.status = censoring.status,
             featurenames = featurenames)


train.obj <- superpc.train(data, type = "survival")
cv.obj <- superpc.cv(train.obj, data)

superpc.plotcv(cv.obj)


## Prefiltering the windows to reduce the dimensionality of the regression
# Keeping 50k windows with the highest frequency, without filtering for the PFS indicator
PLS_data_filter <- studied_patients_data[, .(Patient_ID, PFS_indicator_adjusted)]
PLS_data_filter[, Patient_ID := gsub("NSLC-", "", Patient_ID)]
most_mutated_sw.id <- chrom_frequency[, order(-variant_frequency)][1:50000]
most_mutated_sw.data <- chrom_frequency[most_mutated_sw.id]

pb <- progress_bar$new(format = "[:bar] :percent",
                       total = nrow(most_mutated_sw.data),
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar

for(window in 1:nrow(most_mutated_sw.data)) {
  pb$tick()
  
  if(most_mutated_sw.data[window]$variant_frequency == 0) next
  
  window_name <- paste(c(most_mutated_sw.data[window]$CHROM, most_mutated_sw.data[window]$window_start, most_mutated_sw.data[window]$window_end), collapse=".")
  window_patients <- as.data.table(plyr::count(unlist(str_split(most_mutated_sw.data[window, .(Patients_ID)], ";"))))
  colnames(window_patients) <- c("Patient_ID", window_name)
  if(window %% 1000 == 0) {
    start.time <- Sys.time()
    PLS_data_filter[, glue(window_name) := merge(x=PLS_data_filter[, .(Patient_ID)], y=window_patients, by = "Patient_ID", all.x = TRUE)[[glue(window_name)]]]
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  } else {
    PLS_data_filter[, glue(window_name) := merge(x=PLS_data_filter[, .(Patient_ID)], y=window_patients, by = "Patient_ID", all.x = TRUE)[[glue(window_name)]]]
  }
}


colnames(PLS_data_filter) <- as.character(gsub("\\.", ":", colnames(PLS_data_filter)))

# Event patients are defined as PFS < 2 years with death or relapse
options(scipen=999)
model <- plsr(PFS_indicator_adjusted ~ ., data=PLS_data_filter[, 2:10000], scale=TRUE, validation="CV")
summary(model)
validationplot(model)
validationplot(model, val.type="R2")
















