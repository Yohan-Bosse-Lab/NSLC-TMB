# Incidence-based time-dependent ROC curves
# Reference : https://rstudio-pubs-static.s3.amazonaws.com/3506_36a9509e9d544386bd3e69de30bca608.html
library(risksetROC)
library(tidyverse)
library(ggpubr)
library(data.table)
library(glue)
library(xlsx)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################## Data import and formatting ############################
################################################################################

WGS_and_optimized_regions.data <- fread(file="Results/Feature Selection/studied_patients_data_5k_50k_500k_TMB.csv")

################################################################################
############################## Model building ##################################
################################################################################

# Progression-free survival with pathological stage only
PFS.pstage <- survival::coxph(data = WGS_and_optimized_regions.data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ pathological_stage_refactor)

# Progression-free survival with WGS TMB
PFS.WGS <- survival::coxph(data = WGS_and_optimized_regions.data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low)
PFS.WGS.pstage <- survival::coxph(data = WGS_and_optimized_regions.data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ TMB_high_low + pathological_stage_refactor)

# Progression-free survival with TMB of sliding windows of length 0.005 Mb
PFS.5k <- survival::coxph(data = WGS_and_optimized_regions.data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_5k)
PFS.5k.pstage <- survival::coxph(data = WGS_and_optimized_regions.data, Surv(time=PFS_time, PFS_indicator_adjusted) ~ sw_TMB_cv_5k + pathological_stage_refactor)

################################################################################
################################ ROC curves ####################################
################################################################################

model_list <- list("PFS.pstage", "PFS.WGS", "PFS.WGS.pstage", "PFS.5k", "PFS.5k.pstage")
predict_time_list <- c(0,365,730)

for(model in model_list) {
  for(predict_time in predict_time_list) {
    # Computing ROC curve at time t = 0, 365 or 730 days
    ROC <- with(WGS_and_optimized_regions.data,
                risksetROC(Stime        = PFS_time,
                           status       = PFS_indicator_adjusted,
                           marker       = get(model)$linear.predictor,
                           predict.time = predict_time,
                           plot         = FALSE,
                           method       = "Cox"))
    
    # Extracting AUC data
    AUC <- round(ROC$AUC, digits = 2)
    
    # Editing ROC variable names
    ROC <- data.table("Sensitivity" = ROC$TP, "1 - Specificity" = ROC$FP)
    
    # Saving ROC and AUC data
    assign(glue("{model}.ROC.t{predict_time}"), ROC)
    assign(glue("{model}.AUC.t{predict_time}"), AUC)
  }
}

## Making the ROC curves of models tested above
# t = 0 days
PFS.pstage.AUC.t0
PFS.WGS.pstage.AUC.t0
PFS.5k.pstage.AUC.t0
ROCplot.t0 <- ggplot(mapping = aes(x=`1 - Specificity`, y=`Sensitivity`)) +
  geom_line(data = PFS.5k.pstage.ROC.t0, aes(color = "#E9C46A"), linewidth = 1.2) +
  geom_line(data = PFS.WGS.pstage.ROC.t0, aes(color = "#E63946"), linewidth = 1.2) +
  geom_line(data = PFS.pstage.ROC.t0, aes(color = "#103F59"), linewidth = 1.2) +
  theme_classic() +
  geom_abline(linetype="dashed") +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  ggtitle("t = 0 days") +
  scale_color_manual("",
                     labels = c(glue("Pathological stage, AUC=0.71"),
                                bquote("TMB"["WGS"] ~ " + pathological stage, AUC=0.77"),
                                bquote("TMB"["opt"] ~ "+ pathological stage, AUC=0.94")),
                     values = c("#103F59", "#E63946", "#E9C46A")) +
  labs(x="1 - Specificity", y="Sensitivity") +
  theme(legend.position=c(.62,.12), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 16))

ggsave(ROCplot.t0, file=glue("Results/Model comparison/path_stage-WGS-5k_ROC_0days.png"), width=6, height=6, dpi=300)

# t = 365 days
PFS.pstage.AUC.t365
PFS.WGS.pstage.AUC.t365
PFS.5k.pstage.AUC.t365
ROCplot.t365 <- ggplot(mapping = aes(x=`1 - Specificity`, y=`Sensitivity`)) +
  geom_line(data = PFS.5k.pstage.ROC.t365, aes(color = "#E9C46A"), linewidth = 1.2) +
  geom_line(data = PFS.WGS.pstage.ROC.t365, aes(color = "#E63946"), linewidth = 1.2) +
  geom_line(data = PFS.pstage.ROC.t365, aes(color = "#103F59"), linewidth = 1.2) +
  theme_classic() +
  geom_abline(linetype="dashed") +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  ggtitle("t = 365 days") +
  scale_color_manual("",
                     labels = c(glue("Pathological stage, AUC=0.71"),
                                bquote("TMB"["WGS"] ~ " + pathological stage, AUC=0.77"),
                                bquote("TMB"["opt"] ~ "+ pathological stage, AUC=0.95")),
                     values = c("#103F59", "#E63946", "#E9C46A")) +
  labs(x="1 - Specificity", y="Sensitivity") +
  theme(legend.position=c(.62,.12), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 16))

ggsave(ROCplot.t365, file=glue("Results/Model comparison/path_stage-WGS-5k_ROC_365days.png"), width=6, height=6, dpi=300)

# t = 730 days
PFS.pstage.AUC.t730
PFS.WGS.pstage.AUC.t730
PFS.5k.pstage.AUC.t730
ROCplot.t730 <- ggplot(mapping = aes(x=`1 - Specificity`, y=`Sensitivity`)) +
  geom_line(data = PFS.5k.pstage.ROC.t730, aes(color = "#E9C46A"), linewidth = 1.2) +
  geom_line(data = PFS.WGS.pstage.ROC.t730, aes(color = "#E63946"), linewidth = 1.2) +
  geom_line(data = PFS.pstage.ROC.t730, aes(color = "#103F59"), linewidth = 1.2) +
  theme_classic() +
  geom_abline(linetype="dashed") +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  ggtitle("t = 730 days") +
  scale_color_manual("",
                     labels = c(glue("Pathological stage, AUC=0.7"),
                                bquote("TMB"["WGS"] ~ " + pathological stage, AUC=0.76"),
                                bquote("TMB"["opt"] ~ "+ pathological stage, AUC=0.69")),
                     values = c("#103F59", "#E63946", "#E9C46A")) +
  labs(x="1 - Specificity", y="Sensitivity") +
  theme(legend.position=c(.62,.12), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 16))

ggsave(ROCplot.t730, file=glue("Results/Model comparison/path_stage-WGS-5k_ROC_730days.png"), width=6, height=6, dpi=300)

ROCplot.all <- ggarrange(ROCplot.t0,ROCplot.t365,ROCplot.t730, ncol = 3, nrow = 1, labels = "AUTO")
ggsave(ROCplot.all, file=glue("Results/Model comparison/path_stage-WGS-5k_ROC_all.png"), width=18, height=6, dpi=300)

################################################################################
################################# Archives #####################################
################################################################################

# Other way of showing the progressive AUC score through time before 2 years (730 days)
PFS.pstage.AUC <- risksetAUC(Stime = WGS_and_optimized_regions.data$PFS_time, 
                             status = WGS_and_optimized_regions.data$PFS_indicator_adjusted, 
                             marker = PFS.pstage$linear.predictors,
                             method = "Cox", tmax=730, col="#103F59", lty = 1)

PFS.WGS.pstage.AUC <- risksetAUC(Stime = WGS_and_optimized_regions.data$PFS_time, 
                                 status = WGS_and_optimized_regions.data$PFS_indicator_adjusted, 
                                 marker = PFS.WGS.pstage$linear.predictors,
                                 method = "Cox", tmax=730, plot = F)

PFS.5k.pstage.AUC <- risksetAUC(Stime = WGS_and_optimized_regions.data$PFS_time, 
                                status = WGS_and_optimized_regions.data$PFS_indicator_adjusted, 
                                marker = PFS.5k.pstage$linear.predictors,
                                method = "Cox", tmax=730, col="#E9C46A", plot = F)

lines(PFS.WGS.pstage.AUC$utimes, PFS.WGS.pstage.AUC$AUC, lty=1, col="#E63946")
lines(PFS.5k.pstage.AUC$utimes, PFS.5k.pstage.AUC$AUC, lty=1, col="#E9C46A")



