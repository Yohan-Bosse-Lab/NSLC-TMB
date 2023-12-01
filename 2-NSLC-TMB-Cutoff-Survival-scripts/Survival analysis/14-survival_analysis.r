library(readxl)
library(survival)
library(survminer)
library(lubridate)
library(dplyr)
library(grid)
library(xlsx)
library(lmtest)
library(splines)
library(htmltools)
library(webshot)
library(formattable)
library(AICcmodavg)
library(tidyr)

# *** THIS SCRIPT NEEDS TO BE UPDATED WITH NEW TMB VALUES (AFTER REMOVAL OF dbSNP) ***


export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

# Surival analysis reference: http://www.sthda.com/english/wiki/survival-analysis-basics
# Cox model: http://www.sthda.com/english/wiki/cox-proportional-hazards-model

### Reading combined data table (containing patient clinical characteristics, TMB scores and mutation counts)
NSLC_TMB <- read_xlsx("../../NCI_NeverSmoker_n131_20210812_TMB_124.xlsx") # Dating from 2021

# Removing carcinoid tumors from analyses
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$histology == 4,]
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$Patient_ID == "NSLC-0124",]
NSLC_TMB <- NSLC_TMB[!is.na(NSLC_TMB$normal_sample_ID),] # Removing NA row
nrow(NSLC_TMB)

### Histology and sex formating. Remove "#" comments marks to execute, eases visualisation.
#NSLC_TMB$histology[NSLC_TMB$histology == 2] <- "squamous cell carcinoma (2)"
#NSLC_TMB$histology[NSLC_TMB$histology == 3] <- "adenocarcinoma (3)"
#NSLC_TMB$histology[NSLC_TMB$histology == 4] <- "carcinoid tumor (4)"
#NSLC_TMB$histology[NSLC_TMB$histology == 6] <- "adenosquamous carcinoma (6)"
#NSLC_TMB$histology[NSLC_TMB$histology == 7] <- "sarcomatoid carcinoma (7)"
#NSLC_TMB$sex[NSLC_TMB$sex == 1] <- "male"
#NSLC_TMB$sex[NSLC_TMB$sex == 2] <- "female"
NSLC_TMB["pathological_stage"] <- ifelse(NSLC_TMB$path.stage %in% c("1A1","1A2","1A3","1B"), 1, ifelse(NSLC_TMB$path.stage == "2B", 2, 3))

### The following cutoffs were calculated in the previous script "13-TMB_cutoffs.R".
# WGS TMB median cutoff is set at 0.93 mutation/Mb. High TMB is listed as 1 and low TMB is listed as 0.
NSLC_TMB["WGS_TMB_class_med"] <- ifelse(NSLC_TMB$complete_WGS_TMB > 0.93, "high TMB", "low TMB")
# WGS TMB Youden's cutoff is set at 1.2 mutation/Mb. High TMB is listed as 1 and low TMB is listed as 0.
NSLC_TMB["WGS_TMB_class_you"] <- ifelse(NSLC_TMB$complete_WGS_TMB > 1.2, "high TMB", "low TMB")
# WGS TMB Kappa's cutoff is set at 1.35 mutation/Mb. High TMB is listed as 1 and low TMB is listed as 0.
NSLC_TMB["WGS_TMB_class_kappa"] <- ifelse(NSLC_TMB$complete_WGS_TMB > 1.35, "high TMB", "low TMB")
# WGS TMB Chi-square's cutoff is set at 1.30 mutation/Mb. High TMB is listed as 1 and low TMB is listed as 0.
NSLC_TMB["WGS_TMB_class_chi2"] <- ifelse(NSLC_TMB$complete_WGS_TMB > 1.30, "high TMB", "low TMB")
# WGS TMB generalized cutoff (combination of median and Youden's) is set at 1 mutation/Mb. High TMB is listed as 1 and low TMB is listed as 0.
NSLC_TMB["WGS_TMB_class_gen"] <- ifelse(NSLC_TMB$complete_WGS_TMB > 1, "high TMB", "low TMB")
# WGS TMB dummy cutoff is set at 1.25 mutation/Mb. High TMB is listed as 1 and low TMB is listed as 0. 
# This is to test whether a cutoff between Youden's and Chi-square's would be better.
NSLC_TMB["WGS_TMB_class_dummy"] <- ifelse(NSLC_TMB$complete_WGS_TMB > 1.25, "high TMB", "low TMB")
# % of never-smokers classified with a high TMB:
nrow(NSLC_TMB[NSLC_TMB$WGS_TMB_class_gen == "high TMB",])/nrow(NSLC_TMB)

## Number of patients in need of treatment according to cutoffs:
# Generalized cutoff
length(which(NSLC_TMB$WGS_TMB_class_gen == "high TMB")) # 73 patients in need of treatment (44%)
# Youden cutoff
length(which(NSLC_TMB$WGS_TMB_class_you == "high TMB")) # 39 patients in need of treatment (31.2%)
# 13% of patients rest between a cutoff of 1 and 1.20 mutations/Mb. This region could be implemented as "intermediary risk".

### Models to compare
## Series 1: comparing variables
model1.1 <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + histology, data = NSLC_TMB) # all relevent variables
summary(model1.1)
model1.2 <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age, data = NSLC_TMB) # without histology
summary(model1.2)
# Non linear age reference: https://www.drizopoulos.com/courses/emc/survival%20analysis%20in%20r%20companion
model1.3 <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + ns(age, 3), data = NSLC_TMB) # without histology and with nonlinear age
summary(model1.3)
model1.4 <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + ns(age, 3) + histology, data = NSLC_TMB) # with histology and with nonlinear age
summary(model1.4)
model1.5 <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage, data = NSLC_TMB) # only pathological stage
summary(model1.5)

## Series 2: comparing cutoff models
model2.gen <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_gen), data = NSLC_TMB) # WGS TMB cutoff at 1
summary(model2.gen)
model2.med <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_med), data = NSLC_TMB) # WGS TMB cutoff at 0.93
summary(model2.med)
model2.you <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_you), data = NSLC_TMB) # WGS TMB cutoff at 1.20
summary(model2.you)
model2.kappa <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_kappa), data = NSLC_TMB) # WGS TMB cutoff at 1.35
summary(model2.kappa)
model2.chi2 <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_chi2), data = NSLC_TMB) # WGS TMB cutoff at 1.30
summary(model2.chi2)

### Likelihood ratio test:
# Reference : https://api.rpubs.com/tomanderson_34/lrt
## Series 1: comparing variables
variables_lrtests <- list(
lrtest(model1.1, model1.2), # Inclusion or exclusion of histology does not show substantial accuracy difference (non significant p-value)
lrtest(model1.1, model1.3), 
lrtest(model1.1, model1.4), 
lrtest(model1.1, model1.5), # All variables vs only pathological stage do not show substantial accuracy difference (non significant p-value)
lrtest(model1.2, model1.3), # Linear and non-linear age does not show substantial accuracy difference (non significant p-value)
lrtest(model1.2, model1.4),
lrtest(model1.2, model1.5),
lrtest(model1.3, model1.4),
lrtest(model1.3, model1.5),
lrtest(model1.4, model1.5))

## Series 1 table representation
i = 1
var_pvalue_list <- c()
while (i<length(variables_lrtests)+1) {
  var_pvalue_list <- c(var_pvalue_list, round(variables_lrtests[[i]]$Pr[2], digits = 3))
  i <- i+1
} # Extraction of all p-values from series 1 lrtests

lr.s1_df <- data.frame(Model=c("1", "2", "3", "4", "5"),
                       Lr_test=c(-142.73, -142.75, -141.71, -141.66, -144.58))
colnames(lr.s1_df) <- c(colnames(lr.s1_df)[1], "LogLikelihood")
variables_table <- formattable(lr.s1_df, align=c("l", "c"))
export_formattable(variables_table, "./variable_likelihood_table.png", width = "25%")

# Build a null matrix with the final dimensions (5 rows, 4 columns)
lr.s1_mx <- matrix(0, 5, 4)
# Transform the last matrix into a lower triangular matrix without the diagonal and keep track of the 
lwr.s1_mx <- which(lower.tri(lr.s1_mx, diag = FALSE), arr.ind=TRUE)
# Apply the p-values to the corresponding lower triangular matrix position
lr.s1_mx[lwr.s1_mx] <- as.numeric(var_pvalue_list)
# Change the columns and rows name
rownames(lr.s1_mx) <- lr.s1_df$Model
colnames(lr.s1_mx) <- lr.s1_df$Model[1:4]
# Final p-value of lrtest comparison matrix
t(lr.s1_mx)

# In conclusion, there is no significant accuracy differences between all models. 
# However, model 3 (pathological stage, sex and non-linear age) seems to have the best accuracy.
# ** The final chosen model is the #2, after consultation with PR. It is the usual model used in biomedical 
# research and it seems a good choice since none of the models have any statistical difference between them.

## Series 2: comparing cutoff models
lrtest(model2.gen, model2.med) # Generalized cutoff is significantly more accurate than median cutoff (significant p-value)
lrtest(model2.gen, model2.you) # Generalized cutoff is significantly less accurate than Youden cutoff (significant p-value)
lrtest(model2.med, model2.you) # Median cutoff is significantly less accurate than Youden's cutoff (significant p-value)
lrtest(model2.gen, model2.kappa) # Generalized cutoff is significantly less accurate than Kappa's cutoff (significant p-value)
lrtest(model2.gen, model2.chi2) # Generalized cutoff is significantly less accurate than Chi-square's cutoff (significant p-value)
lrtest(model2.you, model2.kappa) # Youden's cutoff is significantly more accurate than Kappa's cutoff (significant p-value)
lrtest(model2.you, model2.chi2) # Youden's cutoff is significantly more accurate than Chi-square's cutoff (significant p-value)
lrtest(model2.kappa, model2.chi2) # Kappa's cutoff is significantly less accurate than Chi-square's cutoff (significant p-value)

model2.dummy <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_dummy), data = NSLC_TMB) # WGS TMB cutoff at 1.25
summary(model2.dummy)
lrtest(model2.you, model2.dummy) # Youden's cutoff is significantly less accurate than the dummy cutoff (significant p-value)

# Series 2 table representation
lr.s2_df <- data.frame(Statistic=c("Generalized","Median", "Youden", "Kappa", "Chi-square"), 
                    Cutoff=c(1, 0.93, 1.20, 1.35, 1.30),
                    Lr_test=c(-116.62, -117.92, -113.59, -114.21, -113.95))
colnames(lr.s2_df) <- c(colnames(lr.s2_df)[1], "Cutoff (rounded)", "LogLikelihood")
cutoff_table <- formattable(lr.s2_df, align=c("l", "c", "c"))
export_formattable(cutoff_table, "./cutoff_likelihood_table.png", width = "50%")

# Unlike Series 1 lrtests, Series 2 lrtest p-value comparison are all significant at < 2.2e-16 (lrtest function digit limit).



### AIC (Akaike information criterion)
# Reference: https://www.scribbr.com/statistics/akaike-information-criterion/
## Series 1 models comparison
models_series_1 <- list(model1.1, model1.2, model1.3, model1.4, model1.5)
models_series_1.names <- c('path stage + sex + age + histology', 'path stage + sex + age', 'path stage + sex + NL age', 'path stage + sex + NL age + histology', 'path stage only')
aictab(cand.set = models_series_1, modnames = models_series_1.names)
# According to the AIC, model 5 would be the best-fit. However, that model only takes into account one variable.
# Thus, model 2 will be chosen, since it takes into account more variables and is generally the one used in biomed research.

## New model based from model to compare the interaction of variables
model1.2.interaction <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage * sex * age, data = NSLC_TMB)

models_series_AIC <- list(model1.2, model1.2.interaction)
models_series_AIC.names <- c('2', '2 interaction')
aictab(cand.set = models_series_AIC, modnames = models_series_AIC.names)
# The interaction variant of model 2 is not better than the combination variant of model 2 (basic model 2)

## Series 2 models comparison
models_series_2 <- list(model2.med, model2.you, model2.gen, model2.kappa, model2.chi2)
models_series_2.names <- c('median', 'youden', 'generalized', 'kappa', 'chi2')
aictab(cand.set = models_series_2, modnames = models_series_2.names)



### Cox Proportional-Hazards Model (simple representation, see SAS script "TMB_Surv_analysis.sas" for complete Cox plots). 
# The model is stratified by the TMB value of WGS TMB (high or low). The generalized cutoff of 1 mutation/Mb is used (WGS_TMB_class_gen).
cox.multivar <- coxph(Surv(time=survival_months, event=death) ~ pathological_stage + sex + age + strata(WGS_TMB_class_gen), data = NSLC_TMB)
#summary(cox.multivar)
ggsurvplot(data=NSLC_TMB, survfit(cox.multivar), palette = c("#E76F51", "#2A9D8F"),
           ggtheme = theme_minimal(), surv.median.line = "hv", risk.table = TRUE)



### Cox Proportional-Hazards Model for evaluating stage, EGFR positive and driver positive influence on OS survival
NSLC_drivers <- read_xlsx("../../NCI_NeverSmoker_n131_20210812_TMB_92_drivers.xlsx")
NSLC_drivers["TMB.high.low"] <- ifelse(NSLC_drivers$complete_WGS_TMB >= 1.70, "high", "low")
NSLC_drivers$pathological_stage <- ifelse(NSLC_drivers$pathological_stage %in% c("1A1","1A2","1A3","1B"), 1, ifelse(NSLC_drivers$pathological_stage == "2B", 2, 3))

#Model 1 : OS = TMB
summary(coxph(Surv(time_os, VitalStatus == 1) ~ TMB.high.low, data = NSLC_drivers))
summary(coxph(Surv(time_os, VitalStatus == 1) ~ TMB.high.low, data = NSLC_drivers))$waldtest[3]

#Model 2 : OS = TMB + stage
summary(coxph(Surv(time_os, VitalStatus == 1) ~ TMB.high.low + pathological_stage, data = NSLC_drivers))
summary(coxph(Surv(time_os, VitalStatus == 1) ~ TMB.high.low + pathological_stage, data = NSLC_drivers))$waldtest[3]

#Model 3 : OS = TMB + EGFR
#Model 4 : OS = TMB + stage + EGFR
#Model 5 : OS = TMB + drivers
#Model 6 : OS = TMB + stage + drivers


##################################
############ ARCHIVES ############
##################################

# Comparison of sex survival analysis, one curve for males and another for females
fit <- survfit(Surv(time = survival_months, event = death) ~ sex, data = NSLC)
print (fit)

ggsurv <- ggsurvplot(fit,
                     pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF"),
                     xlab = "Time in months", 
                     break.time.by = 30)

ggsurv$table <- ggsurv$table + 
  theme(plot.title = element_text("Patients at Risk"))


# Not comparison of groups, one curve for all patients
fit <- survfit(Surv(time = survival_months, event = death) ~ 1, data = NSLC)
print (fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))