library(readxl)
library(cutpointr)
library(OptimalCutpoints)
library(dplyr)
library(DescTools)
library(formattable)
library(htmltools)
library(webshot)
library(ggplot2)
library(ggpubr)
library(epiR)

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

sensitivity_specificity <- function(matrix){
  # matrix format:
  #     Present  Absent
  # Present   A | B
  #           --|--
  #  Absent   c | D
  # A = true positive
  # D = true negative
  A <- matrix[1,1]
  B <- matrix[1,2]
  C <- matrix[2,1]
  D <- matrix[2,2]
  
  boo <- as.table(matrix(c(A,B,C,D), nrow = 2, byrow = TRUE))
  rval <- epi.tests(boo, conf.level = 0.95)
  print(rval)
}

# Reading combined data table (containing patient clinical characteristics, TMB scores and mutation counts)
NSLC_TMB <- read_xlsx("../../NCI_NeverSmoker_n131_20210812_TMB_92.xlsx") # Dating from 2021

# Removing carcinoid tumors from analyses
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$histology == 4,]
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$Patient_ID == "NSLC-0124",]
NSLC_TMB <- NSLC_TMB[!is.na(NSLC_TMB$normal_sample_ID),] # Removing NA row
nrow(NSLC_TMB)
row.names(NSLC_TMB) <- NULL

##################################
#####       WGS cutoffs      #####
##################################

median(NSLC_TMB$complete_WGS_TMB)
MedianCI(NSLC_TMB$complete_WGS_TMB)
mean(NSLC_TMB$complete_WGS_TMB)
min(NSLC_TMB$complete_WGS_TMB)
max(NSLC_TMB$complete_WGS_TMB)

## Estimating the NSLC WGS TMB cutoff using Youden's index and validation with bootstrap.
cp_you_boot_max <- cutpointr(data=NSLC_TMB, complete_WGS_TMB, VitalStatus, boot_runs = 1000, boot_stratify = TRUE,
                    method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max)
plot(cp_you_boot_max) # integrated plot
# 95% confidence interval for Youden's cutoff:
boot_ci(cp_you_boot_max, optimal_cutpoint, in_bag = FALSE, alpha = 0.05)
# Standard deviation
sd(unlist(cp_you_boot_max$boot[[1]][,1]))

# After several runs, the average optimal cutoff with Youden's Index is ~1.70 ±0.11 mutations/Mb.


## Other tests
# Kappa's test: # *** NEEDS TO BE UPDATED ***
cp_kappa_boot_max <- cutpointr(data=NSLC_TMB, complete_WGS_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                               method = maximize_boot_metric, metric = cohens_kappa) # using bootstrap
summary(cp_kappa_boot_max)
plot(cp_kappa_boot_max) # integrated plot
# 95% confidence interval for Kappa's cutoff:
boot_ci(cp_kappa_boot_max, optimal_cutpoint, in_bag = FALSE, alpha = 0.05)

# Chi-square test:
cp_chi_boot_max <- cutpointr(data=NSLC_TMB, complete_WGS_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                             method = maximize_boot_metric, metric = p_chisquared) # using bootstrap
summary(cp_chi_boot_max)
plot(cp_chi_boot_max) # integrated plot
# 95% confidence interval for Chi-square's cutoff:
boot_ci(cp_chi_boot_max, optimal_cutpoint, in_bag = FALSE, alpha = 0.05)

# Using the library optimal.cutpoints. Youden, Kappa index and Chi-square (MinPvalue) tests are used.
opt_cutp_2 <- optimal.cutpoints(X="complete_WGS_TMB", status="death", tag.healthy=0, methods=c("Youden", "MaxKappa", "MinPvalue"), data=NSLC_TMB)
#plot.optimal.cutpoints(opt_cutp_2) # integrated plot



## Data representation:
# Cutoff tests: 
# *** NEEDS TO BE UPDATED ***
cutoff_df <- data.frame(Statistic=c("Median","Youden", "Kappa", "Chi-square"), 
           Cutoff=c(1.10, 1.37, 1.43, 1.55),
           CI=c("0.98-1.25","1.12-1.65", "1.20-1.79", "0.71-2.89"),
           Value=c("-","0.299", "0.324", "p-value=0.0004"))
colnames(cutoff_df) <- c(colnames(cutoff_df)[1:2], "CI (95%)", colnames(cutoff_df)[4])
cutoff_table <- formattable(cutoff_df, align=c("l", "c", "c", "c"))
export_formattable(cutoff_table, "./cutoff_stats_table.png", width = "50%")

## **UNCOMPLETED** Plotting the relationship between TMB and a cutoff to examine its linearity
ggplot(NSLC_TMB, aes(x=complete_WGS_TMB, y=death)) +
  geom_point() +
  theme_classic()


##################################
#####       WES cutoffs      #####
##################################

median(NSLC_TMB$complete_WES_TMB)
MedianCI(NSLC_TMB$complete_WES_TMB)
mean(NSLC_TMB$complete_WES_TMB)
min(NSLC_TMB$complete_WES_TMB)
max(NSLC_TMB$complete_WES_TMB)

## Estimating the NSLC WES TMB cutoff using Youden's index and validation with bootstrap.
cp_you_boot_max_WES <- cutpointr(data=NSLC_TMB, complete_WES_TMB, VitalStatus, boot_runs = 1000, boot_stratify = TRUE,
                             method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max_WES)
plot(cp_you_boot_max_WES) # integrated plot
# 95% confidence interval for Youden's cutoff:
boot_ci(cp_you_boot_max_WES, optimal_cutpoint, in_bag = FALSE, alpha = 0.05)
# Standard deviation
sd(unlist(cp_you_boot_max_WES$boot[[1]][,1]))

# After several runs, the average optimal cutoff with Youden's Index is ~1.25 ± 0.09 mutation/Mb.


#######################################
#####      NGS panels cutoff      #####
#######################################

# F1CDx cutoff - defined as the median
F1CDX_coff <- median(NSLC_TMB$F1CDx_TMB)

# Illumina TSO500 cutoff - defined as the median
Illu_coff <- median(NSLC_TMB$Illu500_TMB)

# MSKI 468 cutoff - defined as the median
MSKI_coff <- median(NSLC_TMB$MSKI468_TMB)

# NeoPlus cutoff - defined as the median
NeoPlus_coff <- median(NSLC_TMB$NeoPlus_TMB)

# OncoMine cutoff - defined as the median
OncoMine_coff <- median(NSLC_TMB$OncoMineTMB_TMB)

# QIAseq cutoff - defined as the median
QIAseq_coff <- median(NSLC_TMB$QIAseq_TMB)

##################################
#####      General info      #####
##################################

## % of patients with high WGS TMB (over 1.70 mutations/Mb).
round(nrow(NSLC_TMB[NSLC_TMB$complete_WGS_TMB >= 1.70,])/nrow(NSLC_TMB), digits = 2)

## General concordance of high TMB patients between WGS and NGS panels. 
# All NGS panels are combined for this calculation, meaning the final median results compares the general result of NGS panels to WGS.
WGS_high <- NSLC_TMB[NSLC_TMB$complete_WGS_TMB >= 1.70, ]$Patient_ID
F1CDx_high <- NSLC_TMB[NSLC_TMB$F1CDx_TMB >= F1CDX_coff, ]$Patient_ID
Illu_high <- NSLC_TMB[NSLC_TMB$Illu500_TMB >= Illu_coff, ]$Patient_ID
MSKI_high <- NSLC_TMB[NSLC_TMB$MSKI468_TMB >= MSKI_coff, ]$Patient_ID
NeoPlus_high <- NSLC_TMB[NSLC_TMB$NeoPlus_TMB >= NeoPlus_coff, ]$Patient_ID
OncoMine_high <- NSLC_TMB[NSLC_TMB$OncoMineTMB_TMB >= OncoMine_coff, ]$Patient_ID
QIAseq_high <- NSLC_TMB[NSLC_TMB$QIAseq_TMB >= QIAseq_coff, ]$Patient_ID

length(WGS_high) # Base number for NGS panels comparison
WGSxF1CDx <- length(intersect(WGS_high, F1CDx_high))/length(WGS_high)
WGSxIllu <- length(intersect(WGS_high, Illu_high))/length(WGS_high)
WGSxMSKI <- length(intersect(WGS_high, MSKI_high))/length(WGS_high)
WGSxNeoPlus <- length(intersect(WGS_high, NeoPlus_high))/length(WGS_high)
WGSxOncoMine <- length(intersect(WGS_high, OncoMine_high))/length(WGS_high)
WGSxQIAseq <- length(intersect(WGS_high, QIAseq_high))/length(WGS_high)
NGS_panels_WGS_percent <- c(WGSxF1CDx, WGSxIllu, WGSxMSKI, WGSxNeoPlus, WGSxOncoMine, WGSxQIAseq)

median(NGS_panels_WGS_percent); min(NGS_panels_WGS_percent); max(NGS_panels_WGS_percent)

## WGS and isWES TMB ratio.
mean(NSLC_TMB$complete_WGS_TMB)/mean(NSLC_TMB$complete_WES_TMB)
median(NSLC_TMB$complete_WGS_TMB)/median(NSLC_TMB$complete_WES_TMB)
# WGS TMBs are, on average, 1.36-fold higher than isWES TMBs.

## Comparison of each panel to WGS in a contingency matrix. This is for supplementary table 2.
# F1CDx
F1CDx.mx <-  matrix(c(nrow(NSLC_TMB[NSLC_TMB$F1CDx_TMB >= F1CDX_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                      nrow(NSLC_TMB[NSLC_TMB$F1CDx_TMB < F1CDX_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                      nrow(NSLC_TMB[NSLC_TMB$F1CDx_TMB >= F1CDX_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ]),
                      nrow(NSLC_TMB[NSLC_TMB$F1CDx_TMB < F1CDX_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ])),
                    nrow=2, ncol=2)


# Illu (TSO 500)
Illu.mx <- matrix(c(nrow(NSLC_TMB[NSLC_TMB$Illu500_TMB >= Illu_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                    nrow(NSLC_TMB[NSLC_TMB$Illu500_TMB < Illu_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                    nrow(NSLC_TMB[NSLC_TMB$Illu500_TMB >= Illu_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ]),
                    nrow(NSLC_TMB[NSLC_TMB$Illu500_TMB < Illu_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ])),
                  nrow=2, ncol=2)


# MSK-IMPACT
MSK.mx <- matrix(c(nrow(NSLC_TMB[NSLC_TMB$MSKI468_TMB >= MSKI_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                   nrow(NSLC_TMB[NSLC_TMB$MSKI468_TMB < MSKI_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                   nrow(NSLC_TMB[NSLC_TMB$MSKI468_TMB >= MSKI_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ]),
                   nrow(NSLC_TMB[NSLC_TMB$MSKI468_TMB < MSKI_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ])),
                 nrow=2, ncol=2)

# NEOplus
NeoPlus.mx <- matrix(c(nrow(NSLC_TMB[NSLC_TMB$NeoPlus_TMB >= NeoPlus_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                       nrow(NSLC_TMB[NSLC_TMB$NeoPlus_TMB < NeoPlus_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                       nrow(NSLC_TMB[NSLC_TMB$NeoPlus_TMB >= NeoPlus_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ]),
                       nrow(NSLC_TMB[NSLC_TMB$NeoPlus_TMB < NeoPlus_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ])),
                     nrow=2, ncol=2)



# OncoMine
OncoMine.mx <- matrix(c(nrow(NSLC_TMB[NSLC_TMB$OncoMineTMB_TMB >= OncoMine_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                        nrow(NSLC_TMB[NSLC_TMB$OncoMineTMB_TMB < OncoMine_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                        nrow(NSLC_TMB[NSLC_TMB$OncoMineTMB_TMB >= OncoMine_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ]),
                        nrow(NSLC_TMB[NSLC_TMB$OncoMineTMB_TMB < OncoMine_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ])),
                      nrow=2, ncol=2)


# QIAseq TMB
QIAseq.mx <- matrix(c(nrow(NSLC_TMB[NSLC_TMB$QIAseq_TMB >= QIAseq_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                      nrow(NSLC_TMB[NSLC_TMB$QIAseq_TMB < QIAseq_coff & NSLC_TMB$complete_WGS_TMB >= 1.70, ]),
                      nrow(NSLC_TMB[NSLC_TMB$QIAseq_TMB >= QIAseq_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ]),
                      nrow(NSLC_TMB[NSLC_TMB$QIAseq_TMB < QIAseq_coff & NSLC_TMB$complete_WGS_TMB < 1.70, ])), 
                    nrow=2, ncol=2)


## Specificity and sensitivity of contingency matrices done above compared to WGS TMB.
sensitivity_specificity(F1CDx.mx)
sensitivity_specificity(Illu.mx)
sensitivity_specificity(MSK.mx)
sensitivity_specificity(NeoPlus.mx)
sensitivity_specificity(OncoMine.mx)
sensitivity_specificity(QIAseq.mx)

# Chi-square p-values associated with contingency matrices above.
chisq.test(F1CDx.mx)
chisq.test(Illu.mx)
chisq.test(MSK.mx)
chisq.test(NeoPlus.mx)
chisq.test(OncoMine.mx)
chisq.test(QIAseq.mx)

## Youden index of all gene panels
# F1CDx
cp_you_boot_max_F1CDx <- cutpointr(data=NSLC_TMB, F1CDx_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                                 method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max_F1CDx)
# cutoff of 4.5 mut/Mb

# TSO 500
cp_you_boot_max_Illu500 <- cutpointr(data=NSLC_TMB, Illu500_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                                   method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max_Illu500)
# cutoff of 2 mut/Mb

# MSK-IMPACT
cp_you_boot_max_MSK <- cutpointr(data=NSLC_TMB, MSKI468_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                                     method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max_MSK)
# cutoff of 2.98 mut/Mb

# OncoMine
cp_you_boot_max_OncoMine <- cutpointr(data=NSLC_TMB, OncoMineTMB_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                                 method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max_OncoMine)
# cutoff of 2.87 mut/Mb

# QIAseq TMB
cp_you_boot_max_QIAseq <- cutpointr(data=NSLC_TMB, QIAseq_TMB, death, boot_runs = 1000, boot_stratify = TRUE,
                                 method = maximize_boot_metric, metric = youden) # using bootstrap
summary(cp_you_boot_max_QIAseq)
# cutoff of 5.03 mut/Mb

##################################
#####         Plots          #####
##################################

Palette <- c("#a63a9b","#3D87B8","#2a9d8f","#9c818e","#e9c46a","#dc8f41","#Fe5bac","#F50025")

###### Correlation plot between WGS and WES to determine WES cutoff

# WGS cutoff is set at 1.70 mutations/Mb

ggplot(NSLC_TMB, aes(x=complete_WGS_TMB, y=complete_WES_TMB, label = Patient_ID)) +
  annotate("rect", xmin = 0, xmax = 1.70, ymin = -0.16, ymax = 1.20, fill = Palette[2], alpha = 0.7) +
  annotate("rect", xmin = 0, xmax = 1.70, ymin = 1.20, ymax = 20, fill = Palette[6], alpha = 0.7) +
  annotate("rect", xmin = 1.70, xmax = 20, ymin = 1.20, ymax = 20, fill = Palette[4], alpha = 0.7) +
  annotate("rect", xmin = 1.70, xmax = 20, ymin = -0.16, ymax = 1.20, fill = Palette[5], alpha = 0.7) +
  geom_label(x = 3.5, y = 14, label = "WGS TMB cutoff of 1.70 mut/Mb\nas assessed by Youden index", 
             hjust = "left", fontface = 2, label.r = unit(0, "pt"), color = Palette[1]) + # WGS label
  annotate(geom = "curve", x = 3.5, y = 14, xend = 1.9, yend = 13, 
           curvature = .3, arrow = arrow(length = unit(3, "mm")), size = 0.75) + # WGS arrow
  geom_label(x = 7.8, y = 3, label = "Equivalent isWES TMB cutoff \ncorrelates to 1.20 mut/Mb", 
             hjust = "left", fontface = 2, label.r = unit(0, "pt"), color = Palette[3]) + # WES label
  annotate(geom = "curve", x = 7.8, y = 3, xend = 6.5, yend = 1.4, 
           curvature = .3, arrow = arrow(length = unit(3, "mm")), size = 0.75) + # WES arrow
  geom_smooth(method = lm, se = FALSE, color = Palette[8]) +
  geom_point() +
  geom_vline(xintercept = 1.70, color = Palette[1], size = 1, linetype="dashed") + # WGS cutoff
  geom_hline(yintercept = 1.20, color = Palette[3], size = 1, linetype="dashed") + # WES cutoff; calculted from correlation equation y = -0.16 + 0.76x as shown in plot.
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, label.x = 8.5, label.y = 6, color = Palette[8]) +
  theme_classic() +
  scale_x_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  scale_y_continuous(name = "isWES TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"), plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")

ggsave(file="Plots/WGSxisWES_comparison.png", width=6, height=6, dpi=300)

## WGS x isWES contingency matrix.
# Reference : https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot/59119854.
TClass <- factor(c("WGS>=1.70", "WGS>=1.70", "WGS<1.70", "WGS<1.70"))
PClass <- factor(c("isWES>=1.20", "isWES<1.20", "isWES>=1.20", "isWES<1.20"))
Y <- c(nrow(NSLC_TMB[NSLC_TMB$complete_WGS_TMB >= 1.7 & NSLC_TMB$complete_WES_TMB >= 1.20,]),
       nrow(NSLC_TMB[NSLC_TMB$complete_WGS_TMB >= 1.7 & NSLC_TMB$complete_WES_TMB < 1.20,]),
       nrow(NSLC_TMB[NSLC_TMB$complete_WGS_TMB < 1.7 & NSLC_TMB$complete_WES_TMB >= 1.20,]),
       nrow(NSLC_TMB[NSLC_TMB$complete_WGS_TMB < 1.7 & NSLC_TMB$complete_WES_TMB < 1.20,]))
df <- data.frame(TClass, PClass, Y)

# Sensitivity and specificity of isWES TMB compared to WGS TMB.
isWES.mx <- matrix(Y, nrow=2, ncol=2)
sensitivity_specificity(isWES.mx)

# Chi-square p-value
chisq.test(isWES.mx)

ggplot(df, aes(x = TClass, y = PClass)) +
  #("Number of patients categorized\nas above or below WGS\nand isWES TMB cutoff values") +
  geom_rect(fill = Palette[c(4,5,6,2)], xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.7) +
  geom_text(aes(label = sprintf("%1.0f", Y)), color = "black", size = 6) +
  theme_bw() +
  facet_grid(PClass~TClass, scales = "free", space = "free", shrink = FALSE, as.table = FALSE, switch = "y",
             labeller = as_labeller(c("WGS<1.70" = "WGS \n<1.70 mut/Mb", 
                                      "WGS>=1.70" = "WGS \n\u22651.70 mut/Mb",
                                      "isWES<1.20" = "isWES \n<1.20 mut/Mb", 
                                      "isWES>=1.20" = "isWES \n\u22651.20 mut/Mb"))) +
  theme(legend.position = "none", axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(),
        plot.title = element_text(size=16, face = "bold", hjust=0.5), panel.spacing = unit(-0, "lines"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_rect(colour="black", fill = NA),
        strip.text = element_text(size=16, color = "black", hjust = 0.5), strip.text.y.left = element_text(angle = 0, margin = margin(0,0.5,0,0.5, "cm")),
        strip.text.x = element_text(margin = margin(1.5,0,1.5,0, "cm")))

ggsave(file="Plots/WGSxisWES_comparison_table.png", width=5.7, height=6, dpi=300, scale = 0.9)

##################################
#####         Extra          #####
##################################

## Classification of patients if above or below WGS TMB cutoff (1.70 mut/Mb)
write.xlsx(NSLC_TMB[NSLC_TMB$complete_WGS_TMB >= 1.7,c("Patient_ID","ID","complete_WGS_TMB")], "Classification_TMB_92_patients_traitement_over.xlsx", rowNames = FALSE)
write.xlsx(NSLC_TMB[NSLC_TMB$complete_WGS_TMB < 1.7,c("Patient_ID","ID","complete_WGS_TMB")], "Classification_TMB_92_patients_traitement_under.xlsx", rowNames = FALSE)



### CHIP variants adjustments
NSLC_TMB_CHIP <- NSLC_TMB

CHIP_list_1 <- c("NSLC-0018","NSLC-0028","NSLC-0029","NSLC-0036","NSLC-0049",
                 "NSLC-0064","NSLC-0065","NSLC-0081","NSLC-0085","NSLC-0087",
                 "NSLC-0146","NSLC-0155","NSLC-0156")

CHIP_list_2 <- c("NSLC-0020","NSLC-0027")

## Adding CHIP variants TMB changes
# Patients with only one CHIP variant
for (pt in CHIP_list_1) {
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WES_mutation_count <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WES_mutation_count - 1
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$complete_WES_TMB <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WES_mutation_count / 32.8
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WGS_mutation_count <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WGS_mutation_count - 1
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$complete_WGS_TMB <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WGS_mutation_count / 3000
}

# Patients with two CHIP variants
for (pt in CHIP_list_2) {
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WES_mutation_count <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WES_mutation_count - 2
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$complete_WES_TMB <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WES_mutation_count / 32.8
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WGS_mutation_count <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WGS_mutation_count - 2
  NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$complete_WGS_TMB <- NSLC_TMB_CHIP[NSLC_TMB_CHIP$Patient_ID == pt,]$WGS_mutation_count / 3000
}

median(NSLC_TMB_CHIP$complete_WGS_TMB)
MedianCI(NSLC_TMB_CHIP$complete_WGS_TMB)
mean(NSLC_TMB_CHIP$complete_WGS_TMB)
min(NSLC_TMB_CHIP$complete_WGS_TMB)
max(NSLC_TMB_CHIP$complete_WGS_TMB)

# For CHIP variants:
median(NSLC_TMB_CHIP$complete_WES_TMB)
MedianCI(NSLC_TMB_CHIP$complete_WES_TMB)
mean(NSLC_TMB_CHIP$complete_WES_TMB)
min(NSLC_TMB_CHIP$complete_WES_TMB)
max(NSLC_TMB_CHIP$complete_WES_TMB)

##################################
############ ARCHIVES ############
##################################

cp_you_max <- cutpointr(data=NSLC_TMB, complete_WGS_TMB, death,
                    method = maximize_metric, metric = youden)

summary(cp_you_max)

cp_you_min <- cutpointr(data=NSLC_TMB, complete_WGS_TMB, death,
                        method = minimize_metric, metric = youden)

summary(cp_you_min)


