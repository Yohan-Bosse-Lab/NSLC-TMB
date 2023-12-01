
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(survival)
library(finalfit)
library(AICcmodavg)
library(lmtest)
library(gt)
library(webshot2)
library(CPE)
library(pROC)
library(risksetROC)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################## Data import and formatting ############################
################################################################################

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/n82_NSLC_surv_data.xlsx", sheetIndex = 1))

# Fetching and merging all patients variants
VEP_data_all_patients <- fread("Data/n82_NSLC_variants_data.csv", keepLeadingZeros = TRUE)

# Event of interest: progression event before 2 years post-surgery
clin_data[, PFS_indicator_adjusted := ifelse(PFS_indicator == 1 & pfs_two == "<2", 1, 0)]

# The patients taken 
studied_patients_data <- clin_data[!(PFS_indicator == 0 & pfs_two == "<2"), ]

################################################################################
############################## Model building ##################################
################################################################################

# Progression-free survival with basic clinical variables
prog.clinical <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ age + sex + pathological_stage_refactor)
cox.zph(prog.clinical)
summary(prog.clinical)
prog.pstage <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ pathological_stage_refactor)
cox.zph(prog.pstage)
summary(prog.pstage)

## Progression-free survival with continuous WGS TMB
prog.con <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_per_region_genome_TMB)
cox.zph(prog.con)
prog.con.pstage <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_per_region_genome_TMB + pathological_stage_refactor)
cox.zph(prog.con.pstage)
prog.con.pstage.inter <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_per_region_genome_TMB*pathological_stage_refactor)
cox.zph(prog.con.pstage.inter)
prog.con.clinical <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_per_region_genome_TMB + pathological_stage_refactor + age + sex)
cox.zph(prog.con.clinical)
prog.con.clinical.pstage.inter <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_per_region_genome_TMB*pathological_stage_refactor + age + sex)
cox.zph(prog.con.clinical.pstage.inter)

## Progression-free survival with dichotomous WGS TMB
prog.dicho <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low)
cox.zph(prog.dicho)
prog.dicho.pstage <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low + pathological_stage_refactor)
cox.zph(prog.dicho.pstage)
prog.dicho.pstage.inter <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low*pathological_stage_refactor)
cox.zph(prog.dicho.pstage.inter)
prog.dicho.clinical <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low + pathological_stage_refactor + age + sex)
cox.zph(prog.dicho.clinical)
prog.dicho.clinical.pstage.inter <- survival::coxph(data = studied_patients_data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low*pathological_stage_refactor + age + sex)
cox.zph(prog.dicho.clinical.pstage.inter)

## Reporting some of the above models HR in a convenient way
studied_patients_data %>% 
  finalfit("Surv(PFS_time, PFS_indicator_adjusted)", c("TMB_high_low", "pathological_stage_refactor"))

################################################################################
############################# Continuous TMB models ############################
################################################################################

uni.models.con <- list(prog.pstage, prog.con)
uni.mod.names.con <- c("pstage", "continuous")

multi.models.con <- list(prog.clinical, prog.con.pstage, prog.con.pstage.inter, 
                         prog.con.clinical, prog.con.clinical.pstage.inter)
multi.mod.names.con <- c("clinical", "continuous.pstage", "continuous.pstage.interaction",
                         "continuous.clinical", "continuous.pstage.interaction.clinical")

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
                               uni.C.results.con, uni.C.se.con) %>% reduce(full_join, by="Model")

uni.models.results.con <- uni.models.results.con %>% mutate(Model = case_when(Model == "pstage" ~ "pathological stage",
                                                                              Model == "continuous" ~ "TMB"))

multi.models.results.con <- list(AIC.results.multi.con, multi.wald.p.con,
                                 multi.C.results.con, multi.C.se.con) %>% reduce(full_join, by="Model")

multi.models.results.con <- multi.models.results.con %>% mutate(Model = case_when(Model == "continuous.pstage" ~ "TMB + pathological stage",
                                                                                  Model == "clinical" ~ "pathological stage + age + sex",
                                                                                  Model == "continuous.clinical" ~ "TMB + pathological stage + age + sex",
                                                                                  Model == "continuous.pstage.interaction" ~ "TMB * pathological stage",
                                                                                  Model == "continuous.pstage.interaction.clinical" ~ "TMB * pathological stage + age + sex",))

models.results.con <- list(uni.models.results.con, multi.models.results.con) %>% reduce(full_join)


# Likelihood ratio test between the pathological stage model and the pathological stage + TMB model
LRT.con <- round(lmtest::lrtest(prog.pstage, prog.con.pstage)$`Pr(>Chisq)`[2], digits = 3)

# Models comparison table
models.results.con[c(1,2,4,3,6,5,7), .(Model, p, `C-index`, se, AIC)] %>%
  mutate_at(vars(p) , ~(round(., 3))) %>% mutate_at(vars(AIC, `C-index`, se) , ~(round(., 2))) %>% # rounding numeric values
  gt(rowname_col = "Model") %>% cols_align(align = "center", columns = c(p, `C-index`, se, AIC)) %>% # making the table
  cols_label(p = "Wald P") %>%
  tab_stubhead(label = "Model") %>%
  tab_style(style = cell_text(align = "left", indent = px(20)),
            locations = cells_stub()) %>%
  tab_row_group(rows = 3:7, label = "Multivariate") %>%
  tab_row_group(rows = 1:2, label = "Univariate") %>%
  tab_header(title = "Progression-free survival model comparison", subtitle = "Continous whole-genome TMB, n=73") %>%
  tab_footnote(footnote = glue("pathological stage and pathological stage + TMB models Likelihood ratio test p = {LRT.con}"), 
               locations = cells_stub(rows = c("pathological stage", "TMB + pathological stage"))) %>%
  gtsave(filename = "Results/Model comparison/adjusted_PFS_continuous_TMB.docx")



################################################################################
############################# Dichotomous TMB models ###########################
################################################################################

uni.models.dicho <- list(prog.pstage, prog.dicho)
uni.mod.names.dicho <- c("pstage", "dichotomous")

multi.models.dicho <- list(prog.clinical, prog.dicho.pstage, prog.dicho.pstage.inter, 
                           prog.dicho.clinical, prog.dicho.clinical.pstage.inter)
multi.mod.names.dicho <- c("clinical", "dichotomous.pstage", "dichotomous.pstage.interaction",
                           "dichotomous.clinical", "dichotomous.pstage.interaction.clinical")

# Univariate AIC results
AIC.results.uni.dicho <- as.data.table(aictab(cand.set =  uni.models.dicho, modnames = uni.mod.names.dicho))
colnames(AIC.results.uni.dicho)[c(1,3)] <- c("Model", "AIC")

# Multivariate AIC results
AIC.results.multi.dicho <- as.data.table(aictab(cand.set =  multi.models.dicho, modnames = multi.mod.names.dicho))
colnames(AIC.results.multi.dicho)[c(1,3)] <- c("Model", "AIC")

# Univariate Wald test p-values
uni.wald.p.dicho <- as.data.table(setNames(sapply(uni.models.dicho, function(x) summary(x)$waldtest[3]),
                                           nm=uni.mod.names.dicho), keep.rownames = TRUE)
setnames(uni.wald.p.dicho, c("Model", "p"))

# Multivariate Wald test p-values
multi.wald.p.dicho <- as.data.table(setNames(sapply(multi.models.dicho, function(x) summary(x)$waldtest[3]), 
                                             nm=multi.mod.names.dicho), keep.rownames = TRUE)
setnames(multi.wald.p.dicho, c("Model", "p"))

# Univariate concordance index
# Reference : https://rstudio-pubs-static.s3.amazonaws.com/3506_36a9509e9d544386bd3e69de30bca608.html
uni.C.results.dicho <- as.data.table(setNames(sapply(uni.models.dicho, function(x) phcpe(x, CPE.SE = TRUE)$CPE), nm=uni.mod.names.dicho), keep.rownames = TRUE)
setnames(uni.C.results.dicho, c("Model", "C-index"))

uni.C.se.dicho <- as.data.table(setNames(sapply(uni.models.dicho, function(x) phcpe(x, CPE.SE = TRUE)$CPE.SE), nm=uni.mod.names.dicho), keep.rownames = TRUE)
setnames(uni.C.se.dicho, c("Model", "se"))

# Multivariate concordance index
multi.C.results.dicho <- as.data.table(setNames(sapply(multi.models.dicho, function(x) phcpe(x, CPE.SE = TRUE)$CPE), nm=multi.mod.names.dicho), keep.rownames = TRUE)
setnames(multi.C.results.dicho, c("Model", "C-index"))

multi.C.se.dicho <- as.data.table(setNames(sapply(multi.models.dicho, function(x) phcpe(x, CPE.SE = TRUE)$CPE.SE), nm=multi.mod.names.dicho), keep.rownames = TRUE)
setnames(multi.C.se.dicho, c("Model", "se"))

# Preparing the table output
uni.models.results.dicho <- list(AIC.results.uni.dicho, uni.wald.p.dicho,
                                 uni.C.results.dicho, uni.C.se.dicho) %>% reduce(full_join, by="Model")

uni.models.results.dicho <- uni.models.results.dicho %>% mutate(Model = case_when(Model == "pstage" ~ "pathological stage",
                                                                                  Model == "dichotomous" ~ "TMB"))

multi.models.results.dicho <- list(AIC.results.multi.dicho, multi.wald.p.dicho,
                                   multi.C.results.dicho, multi.C.se.dicho) %>% reduce(full_join, by="Model")

multi.models.results.dicho <- multi.models.results.dicho %>% mutate(Model = case_when(Model == "dichotomous.pstage" ~ "TMB + pathological stage",
                                                                                      Model == "clinical" ~ "pathological stage + age + sex",
                                                                                      Model == "dichotomous.clinical" ~ "TMB + pathological stage + age + sex",
                                                                                      Model == "dichotomous.pstage.interaction" ~ "TMB * pathological stage",
                                                                                      Model == "dichotomous.pstage.interaction.clinical" ~ "TMB * pathological stage + age + sex",))

models.results.dicho <- list(uni.models.results.dicho, multi.models.results.dicho) %>% reduce(full_join)


# Likelihood ratio test between the pathological stage model and the pathological stage + TMB model
LRT.dicho <- round(lmtest::lrtest(prog.pstage, prog.dicho.pstage)$`Pr(>Chisq)`[2], digits = 3)

# Models comparison table
models.results.dicho[c(2,1,6,3,5,4,7), .(Model, p, `C-index`, se, AIC)] %>%
  mutate_at(vars(p) , ~(round(., 3))) %>% mutate_at(vars(AIC, `C-index`, se) , ~(round(., 2))) %>% # rounding numeric values
  gt(rowname_col = "Model") %>% cols_align(align = "center", columns = c(p, `C-index`, se, AIC)) %>% # making the table
  cols_label(p = "Wald P") %>%
  tab_stubhead(label = "Model") %>%
  tab_style(style = cell_text(align = "left", indent = px(20)),
            locations = cells_stub()) %>%
  tab_row_group(rows = 3:7, label = "Multivariate") %>%
  tab_row_group(rows = 1:2, label = "Univariate") %>%
  tab_header(title = "Progression-free survival model comparison", subtitle = "Dichotomous whole-genome TMB (< & \u2265 1.70 mut/Mb), n=73") %>%
  tab_footnote(footnote = glue("pathological stage and pathological stage + TMB models Likelihood ratio test p = {LRT.dicho}"), 
               locations = cells_stub(rows = c("pathological stage", "TMB + pathological stage"))) %>%
  gtsave(filename = "Results/Model comparison/adjusted_PFS_dichotomous_TMB.docx")


# Incidence-based time-dependent ROC curves
# Reference : https://rstudio-pubs-static.s3.amazonaws.com/3506_36a9509e9d544386bd3e69de30bca608.html
library(risksetROC)

# ROC curve for dichotomous TMB + pathological stage PFS at 2 years
prog.dicho.pstage.roc.2_years <- with(studied_patients_data,
                                      risksetROC(Stime        = PFS_time,
                                                 status       = PFS_indicator_adjusted,
                                                 marker       = prog.dicho.pstage$linear.predictor,
                                                 predict.time = 730,
                                                 plot         = FALSE,
                                                 method       = "Cox"))
prog.dicho.pstage.AUC.2_years <- round(prog.dicho.pstage.roc.2_years$AUC, digits = 2)
prog.dicho.pstage.roc.2_years <- data.table("Sensitivity" = prog.dicho.pstage.roc.2_years$TP, "1 - Specificity" = prog.dicho.pstage.roc.2_years$FP)

# ROC curve for pathological stage PFS at 2 years
prog.pstage.roc.2_years <- with(studied_patients_data,
                                risksetROC(Stime        = PFS_time,
                                           status       = PFS_indicator_adjusted,
                                           marker       = prog.pstage$linear.predictor,
                                           predict.time = 730,
                                           plot         = FALSE,
                                           method       = "Cox"))
prog.pstage.AUC.2_years <- round(prog.pstage.roc.2_years$AUC, digits = 2)
prog.pstage.roc.2_years <- data.table("Sensitivity" = prog.pstage.roc.2_years$TP, "1 - Specificity" = prog.pstage.roc.2_years$FP)

# Making the ROC curves of both models tested above
ggplot(data = prog.dicho.pstage.roc.2_years, aes(x=`1 - Specificity`, y=`Sensitivity`)) +
  geom_line(aes(color = "#457B9D"), linewidth = 1.2) +
  geom_line(data = prog.pstage.roc.2_years, aes(color = "#E63946"), linewidth = 1.2) +
  theme_classic() +
  geom_abline(linetype="dashed") +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  scale_color_manual("Modèle d'événement de progresssion", 
                     labels = c(glue("TMB + stade pathologique, AUC = {prog.dicho.pstage.AUC.2_years}"), 
                                glue("Stade pathologique, AUC = {prog.pstage.AUC.2_years}")),
                     values = c("#103F59", "#E63946")) +
  labs(x="1 - Spécificité", y="Sensibilité") +
  theme(legend.position=c(.65,.15), 
        legend.title=element_text(size=13), 
        legend.text=element_text(size=11),
        axis.title = element_text(size=13))

ggsave(file=glue("Results/Model comparison/dichotomous_TMB_and_path_stage_ROC_730days.png"), width=5, height=5, dpi=300)


# ROC curve for continuous TMB + pathological stage PFS at 2 years
prog.con.pstage.roc.2_years <- with(studied_patients_data,
                                    risksetROC(Stime        = PFS_time,
                                               status       = PFS_indicator_adjusted,
                                               marker       = prog.con.pstage$linear.predictor,
                                               predict.time = 730,
                                               plot         = FALSE,
                                               method       = "Cox"))
prog.con.pstage.AUC.2_years <- round(prog.con.pstage.roc.2_years$AUC, digits = 2)
prog.con.pstage.roc.2_years <- data.table("Sensitivity" = prog.con.pstage.roc.2_years$TP, "1 - Specificity" = prog.con.pstage.roc.2_years$FP)

# ROC curve for pathological stage PFS at 2 years
prog.pstage.roc.2_years <- with(studied_patients_data,
                                risksetROC(Stime        = PFS_time,
                                           status       = PFS_indicator_adjusted,
                                           marker       = prog.pstage$linear.predictor,
                                           predict.time = 730,
                                           plot         = FALSE,
                                           method       = "Cox"))
prog.pstage.AUC.2_years <- round(prog.pstage.roc.2_years$AUC, digits = 2)
prog.pstage.roc.2_years <- data.table("Sensitivity" = prog.pstage.roc.2_years$TP, "1 - Specificity" = prog.pstage.roc.2_years$FP)

# Making the ROC curves of both models tested above
ggplot(data = prog.con.pstage.roc.2_years, aes(x=`1 - Specificity`, y=`Sensitivity`)) +
  geom_line(aes(color = "#457B9D"), linewidth = 1.2) +
  geom_line(data = prog.pstage.roc.2_years, aes(color = "#E63946"), linewidth = 1.2) +
  theme_classic() +
  geom_abline(linetype="dashed") +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  scale_color_manual("Mod?le d'?v?nement de progresssion", 
                     labels = c(glue("TMB + stade pathologique, AUC = {prog.con.pstage.AUC.2_years}"), 
                                glue("Stade pathologique, AUC = {prog.pstage.AUC.2_years}")),
                     values = c("#457B9D", "#E63946")) +
  labs(x="1 - Sp?cificit?", y="Sensibilit?") +
  theme(legend.position=c(.7,.1))

ggsave(file=glue("Results/Model comparison/continuous_TMB_and_path_stage_ROC_730days.png"), width=5, height=5, dpi=300)


################################################################################
############################# Other statistical tests ##########################
################################################################################

lmtest::lrtest(prog.pstage, prog.dicho)
lmtest::lrtest(prog.pstage, prog.dicho.pstage)
lmtest::lrtest(prog.dicho, prog.dicho.pstage)

