
## TODO
## - Make a list of cancer protectors, cancer drivers and DNA damage repair genes to find mutations
## - For somatic and germline WGS, annotate these genes
## - Make Cox models for drivers with high and low classification
## - Check the count of unique mutations throughout the cohort. Then, separate the cohort into high and low OS to check the ratio of mutations (Sebastien's idea).

library(xlsx)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/n82_NSLC_surv_data.xlsx", sheetIndex = 1))

################################################################################
############################## Overall survival ################################
################################################################################

## TMB as a continuous variable
surv.wg <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB)
surv.wg.pstage <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + pathological_stage_refactor)
surv.exons <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.exons_TMB)
surv.introns <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.introns_TMB)
surv.reg <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.regulatory_TMB)
surv.intergenic <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.intergenic_TMB)

## TMB as a dichotomous variable
surv.di <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_high_low_old)
surv.di.pstage <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_high_low_old + pathological_stage_refactor)

### Verifying Cox P-H model assumptions
# Checking for covariates linearity
plot(predict(surv.wg),
     residuals(surv.wg, type="deviance"))
abline(h=0)
lines(smooth.spline(predict(surv.wg), 
                    residuals(surv.wg, type="deviance")), col = "red")

# Checking for proportional hazards assumption
plot(cox.zph(surv.wg))
abline(h=0, col="red")


### Plotting the survival models
## TMB as a continuous variable
# Infos for plot
surv.wg.p <- round(summary(surv.wg)$wald["pvalue"], digits = 4)
surv.wg.pstage.p <- round(summary(surv.wg.pstage)$wald["pvalue"], digits = 4)
surv.wg.pstage.HR <- round(exp(summary(surv.wg.pstage)$coefficient[1]), digits = 2) # Hazard ratio with Low TMB as reference group
#Plot
surv.wg.plot <- ggsurvplot(survfit(surv.wg, data = surv.dt), 
                                      risk.table = TRUE,
                                      legend.title = "Whole-genome TMB",
                                      ylab = "Survival probability",
                                      risk.table.pos = "in",
                                      conf.int = FALSE)
surv.wg.plot$plot <- surv.wg.plot$plot +
  annotate("text", x=0, y=0.4, hjust = 0,
           label = glue("Wald P-value: {surv.wg.p} \nAdjusted Wald P-value: {surv.wg.pstage.p}\nHazard ratio: {surv.wg.pstage.HR}"))
surv.wg.plot

ggsave("Results/Survival/continuous_OS.png", width=7, height=7, dpi=300)

## TMB as a dichotomous variable
# Infos for plot
surv.di.p <- round(summary(surv.di)$wald["pvalue"], digits = 4)
surv.di.viz <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ strata(TMB_high_low_old)) # for visualization
surv.di.pstage.p <- round(summary(surv.di.pstage)$wald["pvalue"], digits = 4)
surv.di.pstage.HR <- round(exp(-summary(surv.di.pstage)$coefficient[1]), digits = 2) # Hazard ratio with Low TMB as reference group
# Plot
surv.di.plot <- ggsurvplot(survfit(surv.di.viz, data = surv.dt), 
                                      risk.table = TRUE,
                                      legend.title = "Whole-genome TMB",
                                      ylab = "Survival probability",
                                      risk.table.pos = "in")
surv.di.plot$plot <- surv.di.plot$plot +
  annotate("text", x=0, y=0.4, hjust = 0,
           label = glue("Wald P-value: {surv.di.p} \nAdjusted Wald P-value: {surv.di.pstage.p}\nHazard ratio: {surv.di.pstage.HR}"))
surv.di.plot

ggsave("Results/Survival/dicho_OS.png", width=7, height=7, dpi=300)

################################################################################
############################## Recurrence free #################################
################################################################################

## TMB as a continuous variable
surv.wg.reccurence <- coxph(data = surv.dt, Surv(time=time_RpFS, RpFS_indicator) ~ TMB_per_region.genome_TMB)
surv.wg.reccurence.pstage <- coxph(data = surv.dt, Surv(time=time_RpFS, RpFS_indicator) ~ TMB_per_region.genome_TMB + pathological_stage_refactor)

## TMB as a dichotomous variable
surv.di.reccurence <- coxph(data = surv.dt, Surv(time=time_RpFS, RpFS_indicator) ~ TMB_high_low_old)
surv.di.reccurence.pstage <- coxph(data = surv.dt, Surv(time=time_RpFS, RpFS_indicator) ~ TMB_high_low_old + pathological_stage_refactor)


### Plotting the survival models
## TMB as a continuous variable
# Infos for plot
surv.wg.reccurence.p <- round(summary(surv.wg.reccurence)$wald["pvalue"], digits = 4) # Wald p-value
surv.wg.reccurence.pstage.p <- round(summary(surv.wg.reccurence.pstage)$wald["pvalue"], digits = 4) # Wald p-value
surv.wg.reccurence.pstage.HR <- round(exp(summary(surv.wg.reccurence.pstage)$coefficient[1]), digits = 2) # Hazard ratio with Low TMB as reference group
# Plot
surv.wg.reccurence.plot <- ggsurvplot(survfit(surv.wg.reccurence, data = surv.dt), 
                                      risk.table = TRUE,
                                      legend.title = "Whole-genome TMB",
                                      ylab = "Relapse-free probability",
                                      risk.table.pos = "in",
                                      conf.int = FALSE)
surv.wg.reccurence.plot$plot <- surv.wg.reccurence.plot$plot +
  annotate("text", x=0, y=0.4, hjust = 0,
           label = glue("Wald P-value: {surv.wg.reccurence.p} \nAdjusted Wald P-value: {surv.wg.reccurence.pstage.p}\nHazard ratio: {surv.wg.reccurence.pstage.HR}"))
surv.wg.reccurence.plot

ggsave("Results/Survival/continuous_RpFS.png", width=7, height=7, dpi=300)

## TMB as a dichotomous variable
# Infos for plot
surv.di.reccurence.p <- round(summary(surv.di.reccurence)$wald["pvalue"], digits = 4) # Wald p-value
surv.di.reccurence.viz <- coxph(data = surv.dt, Surv(time=time_RpFS, RpFS_indicator) ~ strata(TMB_high_low_old)) # stratified for visualization
surv.di.reccurence.pstage.p <- round(summary(surv.di.reccurence.pstage)$wald["pvalue"], digits = 4) # Wald p-value
surv.di.reccurence.pstage.HR <- round(exp(-summary(surv.di.reccurence.pstage)$coefficient[1]), digits = 2) # Hazard ratio with Low TMB as reference group

# Plot
surv.di.reccurence.plot <- ggsurvplot(survfit(surv.di.reccurence.viz, data = surv.dt), 
           risk.table = TRUE,
           legend.title = "Whole-genome TMB",
           ylab = "Relapse-free probability",
           risk.table.pos = "in")
surv.di.reccurence.plot$plot <- surv.di.reccurence.plot$plot +
  annotate("text", x=0, y=0.4, hjust = 0,
           label = glue("Wald P-value: {surv.di.reccurence.p} \nAdjusted Wald P-value: {surv.di.reccurence.pstage.p}\nHazard ratio: {surv.di.reccurence.pstage.HR}"))
surv.di.reccurence.plot

ggsave("Results/Survival/dicho_RpFS.png", width=7, height=7, dpi=300)

################################################################################
########################### Stratified recurrence ##############################
################################################################################

### Survival analysis with stratified recurrence (< and >= 2 years)
## TMB as a continuous variable
# Model
surv.wg.reccurence_two <- coxph(data = surv.dt[!is.na(recurrence.two)], Surv(time=time_RpFS, RpFS_indicator) ~ TMB_per_region.genome_TMB + strata(recurrence.two))

# Plot
surv.wg.reccurence_two.plot <- ggsurvplot(survfit(surv.wg.reccurence_two),
           data = surv.dt[!is.na(recurrence.two)],
           risk.table = TRUE,
           legend.title = "Time to recurrence",
           ylab = "Relapse-free probability",
           risk.table.pos = "in",
           legend.labs = c("\u2265 2 years", "< 2 years"))

surv.wg.reccurence_two.plot$plot <- surv.wg.reccurence_two.plot$plot +
  geom_vline(xintercept = 730)

surv.wg.reccurence_two.plot
ggsave("Results/Survival/continuous_strata_2_RpFS.png", width=7, height=7, dpi=300)


### Survival analysis with stratified recurrence (< and >= 5 years)
## TMB as a continuous variable
# Model
surv.wg.reccurence_five <- coxph(data = surv.dt[!is.na(recurrence.five)], Surv(time=time_RpFS, RpFS_indicator) ~ TMB_per_region.genome_TMB + strata(recurrence.five))

# Plot
surv.wg.reccurence_five.plot <- ggsurvplot(survfit(surv.wg.reccurence_five),
                                          data = surv.dt[!is.na(recurrence.five)],
                                          risk.table = TRUE,
                                          legend.title = "Time to recurrence",
                                          ylab = "Relapse-free probability",
                                          risk.table.pos = "in",
                                          legend.labs = c("\u2265 5 years", "< 5 years"))

surv.wg.reccurence_five.plot$plot <- surv.wg.reccurence_five.plot$plot +
  geom_vline(xintercept = 1825)

surv.wg.reccurence_five.plot
ggsave("Results/Survival/continuous_strata_5_RpFS.png", width=7, height=7, dpi=300)


################################################################################
################################## Archives ####################################
################################################################################
#
# model.cox.simple <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB)
# 
# model.cox.simple.spline <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ pspline(TMB_per_region.genome_TMB, df=4))
# ptemp <- termplot(model.cox.simple.spline, se=TRUE, plot=FALSE)
# 
# model.cox.pstage <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + pathological_stage_refactor)
# 
# model.cox.pstage.strata <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + strata(pathological_stage_refactor))
# 
# model.cox.age <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + age)
# 
# model.cox.pstage_age <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + age + pathological_stage_refactor)
# 
# anova(model.cox.simple, model.cox.pstage, test="LRT") # The two models do not seem different, therefore we can exclude pathological stage
# 
# ## Check for linearity for the simple model
# # Using Martingale residuals
# plot(predict(model.cox.simple), 
#      residuals(model.cox.simple, type="martingale"))
# abline(h=0)
# lines(smooth.spline(predict(model.cox.simple), 
#                     residuals(model.cox.simple, type="martingale")), col = "red")
# 
# # Using deviance residuals
# plot(predict(model.cox.simple), 
#      residuals(model.cox.simple, type="deviance"))
# abline(h=0)
# lines(smooth.spline(predict(model.cox.simple), 
#                     residuals(model.cox.simple, type="deviance")), col = "red")
# 
# 
# 
# ## Check for proportional hazards assumption for the simple model
# plot(cox.zph(model.cox.simple))
# abline(h=0, col="red")
# 
# ## Check for proportional hazards assumption for the model with pathological stage
# plot(cox.zph(model.cox.pstage)[1])
# abline(h=0, col="red")
# plot(cox.zph(model.cox.pstage)[2]) # The hazards are not proportional for the pathological stage
# abline(h=0, col="red")
# 
# ## Check for proportional hazards assumption for the model with stratified pathological stage
# plot(cox.zph(model.cox.pstage.strata))
# abline(h=0, col="red")




