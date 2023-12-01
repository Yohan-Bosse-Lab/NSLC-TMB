
library(tidyverse)
library(data.table)
library(gghalves)
library(glue)
library(xlsx)
library(progress)
library(survival)
library(finalfit)
library(GenomicRanges)
library(AICcmodavg)
library(webshot2)
library(gt)
library(CPE)
library(pROC)
library(risksetROC)

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

# Fetching and merging all patients variants
VEP_data_all_patients <- fread("Data/n82_NSLC_variants_data.csv", keepLeadingZeros = TRUE)


final_sw.boot <- fread(file=glue("Results/Feature Selection/boot_{Mb_length}Mb_final_features_all_patients.csv"), header = FALSE)
final_sw.cv <- fread(file=glue("Results/Feature Selection/cv_{Mb_length}Mb_final_features_all_patients.csv"), header = FALSE)

################################################################################
##################### TMB of selected sliding windows ##########################
################################################################################

## No need to run this section again
# 
# # Bootstrap sliding windows data for TMB
# final_sw.boot[, c("CHROM", "start", "end") := tstrsplit(V1, ":", fixed = TRUE)]
# final_sw.list.boot <- GRanges(final_sw.boot) %>% reduce() %>% as.data.table()
# final_sw.size.boot <- sum(final_sw.list.boot$width) # size in number of base pairs
# 
# # Cross-validation sliding windows data for TMB
# final_sw.cv[, c("CHROM", "start", "end") := tstrsplit(V1, ":", fixed = TRUE)]
# final_sw.list.cv <- GRanges(final_sw.cv) %>% reduce() %>% as.data.table()
# final_sw.size.cv <- sum(final_sw.list.cv$width) # size in number of base pairs
# 
# # Progression bar for visualizing patients being processed
# pb <- progress_bar$new(format = "[:bar] :percent",
#                        total = length(patient_list),
#                        complete = "=",   # Completion bar character
#                        incomplete = "-", # Incomplete bar character
#                        current = ">",    # Current bar character
#                        clear = FALSE,    # If TRUE, clears the bar when finish
#                        width = 100)      # Width of the progress bar
# 
# for(patient in patient_list) {
#   pb$tick()
#   
#   ## Extracting mutations in selected sliding windows
#   # For boostrap
#   sw_mutations.boot <- rbindlist(lapply(1:nrow(final_sw.list.boot), 
#                                         function(x) (VEP_data_all_patients[Patient_ID == patient &
#                                                                              CHROM == final_sw.list.boot[x, seqnames] &
#                                                                              POS >= final_sw.list.boot[x, start] &
#                                                                              END_POS <= final_sw.list.boot[x, end]])))
#   assign(glue("sw_mutations.boot.{Mb_length}.{patient}"), sw_mutations.boot)
#   
#   # For cross-validation
#   sw_mutations.cv <- rbindlist(lapply(1:nrow(final_sw.list.cv), 
#                                       function(x) (VEP_data_all_patients[Patient_ID == patient &
#                                                                            CHROM == final_sw.list.cv[x, seqnames] &
#                                                                            POS >= final_sw.list.cv[x, start] &
#                                                                            END_POS <= final_sw.list.cv[x, end]])))
#   assign(glue("sw_mutations.cv.{Mb_length}.{patient}"), sw_mutations.cv)
#   
#   # Calculating the TMB of the selected sliding windows
#   # For boostrap
#   sw_TMB.boot <- round(nrow(sw_mutations.boot) / (final_sw.size.boot / 1e6), digits = 3)
#   assign(glue("sw_TMB.boot.{Mb_length}.{patient}"), sw_TMB.boot)
#   
#   # For cross-validation
#   sw_TMB.cv <- round(nrow(sw_mutations.cv) / (final_sw.size.cv / 1e6), digits = 3)
#   assign(glue("sw_TMB.cv.{Mb_length}.{patient}"), sw_TMB.cv)
# }
# 
# sw_TMB.boot.5k.list <- unlist(lapply(objects()[grep(glue("sw_TMB.boot.{Mb_length}"), objects())], function(x) get(x)))
# sw_TMB.cv.5k.list <- unlist(lapply(objects()[grep(glue("sw_TMB.cv.{Mb_length}"), objects())], function(x) get(x)))
# 
# # studied_patients_data[, `:=`(sw_TMB_boot_5k = sw_TMB.boot.5k.list, sw_TMB_cv_5k = sw_TMB.cv.5k.list)]
# 
# fwrite(studied_patients_data, file = "Results/Feature Selection/studied_patients_data_5k_50k_500k_TMB.csv")

################################################################################
############################## Model building ##################################
################################################################################

studied_patients_data <- fread("Results/Feature Selection/studied_patients_data_5k_50k_500k_TMB.csv")

# Progression-free survival with basic clinical variables. Previously tested in script 7-82_NSLC_model_comparison.r
prog.clinical <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ age + sex + pathological_stage_refactor)
prog.pstage <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ pathological_stage_refactor)

## Progression-free survival with dichotomous WGS TMB
prog.dicho.pstage <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low + pathological_stage_refactor)

## Progression-free survival with selected sliding windows TMB
# For sliding windows of length 0.005 Mb
prog.sw_TMB.boot.5k <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_boot_5k)
prog.sw_TMB.boot.5k.path <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_boot_5k + pathological_stage_refactor)
prog.sw_TMB.cv.5k <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_5k)
prog.sw_TMB.cv.5k.path <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_5k + pathological_stage_refactor)

# For sliding windows of length 0.05 Mb
prog.sw_TMB.boot.50k <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_boot_50k)
prog.sw_TMB.cv.50k <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_50k)
prog.sw_TMB.cv.50k.path <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_50k + pathological_stage_refactor)

# For sliding windows of length 0.5 Mb
prog.sw_TMB.boot.500k <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_boot_500k)
prog.sw_TMB.cv.500k <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_500k)

## Reporting some of the above models HR in a convenient way
studied_patients_data %>% 
  finalfit("Surv(PFS_time, PFS_indicator_adjusted)", c("sw_TMB_cv_5k", "pathological_stage_refactor"))

################################################################################
############################# Other statistical tests ##########################
################################################################################

lmtest::lrtest(prog.pstage, prog.sw_TMB.cv.5k.path)
lmtest::lrtest(prog.pstage, prog.sw_TMB.cv.50k.path)
lmtest::lrtest(prog.sw_TMB.cv.5k, prog.sw_TMB.cv.5k.path)
lmtest::lrtest(prog.sw_TMB.cv.50k, prog.sw_TMB.cv.50k.path)

################################################################################
################################# Output models ################################
################################################################################

uni.models.con <- list(prog.pstage, prog.sw_TMB.boot.5k, prog.sw_TMB.cv.5k, prog.sw_TMB.cv.50k)
uni.mod.names.con <- c("pstage", "bootstrap.5k", "cv.5k", "cv.50k")

multi.models.con <- list(prog.sw_TMB.boot.5k.path, prog.sw_TMB.cv.5k.path, prog.sw_TMB.cv.50k.path)
multi.mod.names.con <- c("bootstrap.5k.pstage", "cv.5k.pstage", "cv.50k.pstage")

# Univariate AIC results
AIC.results.uni.con <- as.data.table(aictab(cand.set =  uni.models.con, modnames = uni.mod.names.con))
colnames(AIC.results.uni.con)[c(1,3)] <- c("Model", "AIC")

# Multivariate AIC results
AIC.results.multi.con <- as.data.table(aictab(cand.set =  multi.models.con, modnames = multi.mod.names.con))
colnames(AIC.results.multi.con)[c(1,3)] <- c("Model", "AIC")

# Univariate Wald test p-values
uni.wald.p.con <- as.data.table(setNames(sapply(uni.models.con, function(x) summary(x)$waldtest[3]),
                                         nm=uni.mod.names.con), keep.rownames = TRUE)
setnames(uni.wald.p.con, c("Model", "p"))

# Multivariate Wald test p-values
multi.wald.p.con <- as.data.table(setNames(sapply(multi.models.con, function(x) summary(x)$waldtest[3]), 
                                           nm=multi.mod.names.con), keep.rownames = TRUE)
setnames(multi.wald.p.con, c("Model", "p"))

# Univariate concordance index
# Reference : https://rstudio-pubs-static.s3.amazonaws.com/3506_36a9509e9d544386bd3e69de30bca608.html
uni.C.results.con <- as.data.table(setNames(sapply(uni.models.con, function(x) phcpe(x, CPE.SE = TRUE)$CPE), nm=uni.mod.names.con), keep.rownames = TRUE)
setnames(uni.C.results.con, c("Model", "C-index"))

uni.C.se.con <- as.data.table(setNames(sapply(uni.models.con, function(x) phcpe(x, CPE.SE = TRUE)$CPE.SE), nm=uni.mod.names.con), keep.rownames = TRUE)
setnames(uni.C.se.con, c("Model", "se"))

# Multivariate concordance index
multi.C.results.con <- as.data.table(setNames(sapply(multi.models.con, function(x) phcpe(x, CPE.SE = TRUE)$CPE), nm=multi.mod.names.con), keep.rownames = TRUE)
setnames(multi.C.results.con, c("Model", "C-index"))

multi.C.se.con <- as.data.table(setNames(sapply(multi.models.con, function(x) phcpe(x, CPE.SE = TRUE)$CPE.SE), nm=multi.mod.names.con), keep.rownames = TRUE)
setnames(multi.C.se.con, c("Model", "se"))

# Preparing the table output
uni.models.results.con <- list(AIC.results.uni.con, uni.wald.p.con,
                               uni.C.results.con, uni.C.se.con) %>% purrr::reduce(full_join, by="Model")

uni.models.results.con <- uni.models.results.con %>% mutate(Model = case_when(Model == "pstage" ~ "pathological stage",
                                                                              Model == "bootstrap.5k" ~ "0.005Mb sliding window bootstrap TMB",
                                                                              Model == "cv.5k" ~ "0.005Mb sliding window cross-validation TMB",
                                                                              Model == "cv.50k" ~ "0.05Mb sliding window cross-validation TMB"))

multi.models.results.con <- list(AIC.results.multi.con, multi.wald.p.con,
                                 multi.C.results.con, multi.C.se.con) %>% purrr::reduce(full_join, by="Model")

multi.models.results.con <- multi.models.results.con %>% mutate(Model = case_when(Model == "bootstrap.5k.pstage" ~ "0.005Mb sliding window bootstrap TMB + pathological stage",
                                                                                  Model == "cv.5k.pstage" ~ "0.005Mb sliding window cross-validation TMB + pathological stage",
                                                                                  Model == "cv.50k.pstage" ~ "0.05Mb sliding window cross-validation TMB + pathological stage"))

models.results.con <- list(uni.models.results.con, multi.models.results.con) %>% purrr::reduce(full_join)


# Likelihood ratio test between the pathological stage model and the pathological stage + 5k cv TMB model
LRT.con <- lmtest::lrtest(prog.pstage, prog.sw_TMB.cv.5k.path)$`Pr(>Chisq)`[2]

# Models comparison table
models.results.con[c(4, 1:3, 5:7), .(Model, p, `C-index`, se, AIC)] %>%
  mutate_at(vars(p) , ~(round(., 3))) %>% mutate_at(vars(AIC, `C-index`, se) , ~(round(., 2))) %>% # rounding numeric values
  gt(rowname_col = "Model") %>% cols_align(align = "center", columns = c(p, `C-index`, se, AIC)) %>% # making the table
  cols_label(p = "Wald P") %>%
  tab_stubhead(label = "Model") %>%
  tab_style(style = cell_text(align = "left", indent = px(20)),
            locations = cells_stub()) %>%
  tab_row_group(rows = 5:7, label = "Multivariate") %>%
  tab_row_group(rows = 1:4, label = "Univariate") %>%
  tab_header(title = "Progression-free survival model comparison", subtitle = "Region optimized TMB, n=73") %>%
  tab_footnote(footnote = glue("pathological stage and pathological stage + 0.005Mb CV TMB Likelihood ratio test p = {LRT.con}"),
               locations = cells_stub(rows = c("pathological stage", "0.005Mb sliding window cross-validation TMB + pathological stage"))) %>%
  gtsave(filename = "Results/Feature Selection/adjusted_PFS_region-optimized_TMB.docx")


