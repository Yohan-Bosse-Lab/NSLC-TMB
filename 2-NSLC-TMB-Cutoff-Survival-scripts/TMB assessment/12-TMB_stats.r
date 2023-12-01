# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

library(readxl)
library(ggplot2)
library(ggpubr)
library(grid)
library(ggrepel)
library(ggbreak)
library(RColorBrewer)
library(cowplot)
library(blandr)
library(tibble)
library(dplyr)
library(glue)
library(xlsx)
library(forcats)
library(stringr)
library(magrittr)
  
# ** INFO : isWES means in silico WES

# Importing several patient clinical characteristics and TMB score files.
#NSLC <- read_xlsx("../../NCI_NeverSmoker_n131_20181207_share.xlsx", range = "A1:V132") # Dating from 2018
NSLC_TMB <- read_xlsx("../../NCI_NeverSmoker_n131_20210812_TMB_124.xlsx") # Dating from 2021

# Removing carcinoid tumors from analyses
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$histology == 4,]
NSLC_TMB <- NSLC_TMB[!NSLC_TMB$Patient_ID == "NSLC-0124",]
NSLC_TMB <- NSLC_TMB[!is.na(NSLC_TMB$normal_sample_ID),] # Removing NA row
nrow(NSLC_TMB)

# Export list of patients
#NSLC_TMB <- NSLC_TMB %>% transform(Patient_ID=str_replace(Patient_ID,"NSLC-",""))
#write.csv(NSLC_TMB$Patient_ID, "92_patients_list.csv", row.names = FALSE)

#cols.names <- c("classic_WES_TMB", "complete_WES_TMB", "nc_WGS_TMB", "complete_WGS_TMB", "F1CDx_TMB", "Illu500_TMB", "MSKI468_TMB", "NeoPlus_TMB", "OncoMineTMB_TMB", "QIAseq_TMB", "WES_mutation_count", "WGS_mutation_count")
#NSLC_TMB[cols.names] <- sapply(NSLC_TMB[cols.names], as.character)
#NSLC_TMB[cols.names] <- sapply(NSLC_TMB[cols.names], as.numeric)

# Color theme used in the plots
Palette <- c("#a63a9b","#3D87B8","#2a9d8f","#9c818e","#e9c46a","#dc8f41","#Fe5bac","#F50025")

# Function to identify outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Box plot of all TMB measures, with multiple NGS panels
# Formatting of new dataframe for a boxplot graph of TMBs
boxp_df <- reshape2::melt(NSLC_TMB[,c("Patient_ID", "complete_WES_TMB", "complete_WGS_TMB", "F1CDx_TMB", "Illu500_TMB", "MSKI468_TMB", "NeoPlus_TMB", "OncoMineTMB_TMB", "QIAseq_TMB")], id.vars="Patient_ID", variable.name = "TMB_type", value.name = "TMB")
boxp_df$Patient_ID<-as.numeric(gsub("NSLC-","",as.character(boxp_df$Patient_ID)))
boxp_df$TMB_type <- as.character(boxp_df$TMB_type)
boxp_df$TMB_type[boxp_df$TMB_type == "complete_WGS_TMB"] <- "WGS\n(3000 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "complete_WES_TMB"] <- "isWES\n(32.8 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "F1CDx_TMB"] <- "F1CDx\n(324 genes, 0.8 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "Illu500_TMB"] <- "TSO500\n(523 genes, 1.33 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "MSKI468_TMB"] <- "MSK-IMPACT\n(468 genes, 1.14 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "NeoPlus_TMB"] <- "NEOPlus RUO assay\n(340 genes, 1.1 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "OncoMineTMB_TMB"] <- "Oncomine TML\n(409 genes, 1.2 Mb)"
boxp_df$TMB_type[boxp_df$TMB_type == "QIAseq_TMB"] <- "QIAseq\n(486 genes, 1.33 Mb)"
boxp_df$TMB_type <- factor(boxp_df$TMB_type, levels = c("WGS\n(3000 Mb)", "isWES\n(32.8 Mb)", "F1CDx\n(324 genes, 0.8 Mb)", 
                                                        "TSO500\n(523 genes, 1.33 Mb)", "MSK-IMPACT\n(468 genes, 1.14 Mb)", 
                                                        "NEOPlus RUO assay\n(340 genes, 1.1 Mb)", "Oncomine TML\n(409 genes, 1.2 Mb)", 
                                                        "QIAseq\n(486 genes, 1.33 Mb)"))

boxp_df <- boxp_df %>% tibble::rownames_to_column(var="outlier") %>% group_by(TMB_type) %>% mutate(is_outlier=ifelse(is_outlier(TMB), TMB, as.numeric(NA)))
boxp_df$outlier[which(is.na(boxp_df$is_outlier))] <- as.numeric(NA)
boxp_df$outlier[which(!is.na(boxp_df$is_outlier))] <- as.character(boxp_df$Patient_ID[!is.na(boxp_df$is_outlier)])

top <- ggplot(boxp_df[boxp_df$TMB>=7.5,], aes(x=TMB_type, y=TMB, color=TMB_type)) + 
  #ggtitle("Comparison of TMB values between sequencing methods\n(n = 93 patients)") +
  geom_jitter(shape=21, width = 0.25, size = 2) +
  geom_point(shape=2, size=2.5, color="black") +
  geom_text_repel(data = boxp_df[boxp_df$TMB >= 15,], aes(label=outlier), na.rm = T, nudge_x = 0.35, direction = "y", force = 1) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = "dashed"), 
        panel.grid.minor.y = element_line(colour = "black", size = 0.5, linetype = "dotted"),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.line.x = element_blank(),
        plot.margin = unit(c(0.5, 0.1, 0, 0.5), "cm")) +
  scale_color_manual(values = Palette, "Sequencing method") +
  scale_y_continuous(breaks = seq(10,35,5), minor_breaks = seq(5,35,2.5), expand=c(-0.01,0)) + 
  xlab(NULL) + 
  ylab(NULL) + 
  coord_cartesian(ylim=c(15,35.5), clip = "off")

bot <- ggplot(boxp_df, aes(x=TMB_type, y=TMB, color=TMB_type)) + 
  geom_boxplot(notch = FALSE, outlier.shape = 2, outlier.size = 2.5, color = "black") +
  geom_jitter(shape=21, width = 0.25, size = 2) +
  geom_text_repel(data = boxp_df[boxp_df$TMB < 15,], aes(label=outlier), na.rm = T, nudge_x = 0.35, direction = "y", force = 1) +
  stat_summary(fun=mean, geom="point", shape=18, size=5, color = "black") +
  geom_hline(yintercept=0, size = 0.75) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = "dashed"), 
        panel.grid.minor.y = element_line(colour = "black", size = 0.5, linetype = "dotted"),
        plot.title = element_text(hjust = 0.5, size = 16),   
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14, hjust = 0.9),
        legend.position = "none",
        plot.margin = unit(c(-0.1, 0.1, 0.5, 0.5), "cm")) +
  scale_color_manual(values = Palette, "Sequencing method") +
  scale_y_continuous(breaks = seq(0,15,5), minor_breaks = seq(0,15,2.5), expand=c(0.005,0)) + 
  xlab("Sequencing method") + 
  ylab("TMB (mutations/Mb)") + 
  coord_cartesian(ylim=c(0,15))

plot_grid(top, bot, rel_heights=c(0.33, 0.66), ncol = 1, align = "v")
ggsave(file="TMB graphs/TMB_boxplot.png", width=12, height=11, dpi=300)


###### Box plot of all TMB measures, with only FoundationOne CDx and Illumina 500 panels 
# Formatting of new dataframe for a boxplot graph of TMBs
boxp_df.2 <-reshape2::melt(NSLC_TMB[,c("Patient_ID", "complete_WES_TMB", "complete_WGS_TMB", "F1CDx_TMB", "Illu500_TMB")], id.vars="Patient_ID", variable.name = "TMB_type", value.name = "TMB")
boxp_df.2$Patient_ID<-as.numeric(gsub("NSLC-","",as.character(boxp_df.2$Patient_ID)))
boxp_df.2$TMB_type <- as.character(boxp_df.2$TMB_type)
boxp_df.2$TMB_type[boxp_df.2$TMB_type == "complete_WGS_TMB"] <- "WGS\n(3000 Mb)"
boxp_df.2$TMB_type[boxp_df.2$TMB_type == "complete_WES_TMB"] <- "isWES\n(33 Mb)"
boxp_df.2$TMB_type[boxp_df.2$TMB_type == "F1CDx_TMB"] <- "F1CDx assay\n(324 genes, 2.2 Mb)"
boxp_df.2$TMB_type[boxp_df.2$TMB_type == "Illu500_TMB"] <- "TSO500\n(523 genes, 1.94 Mb)"
boxp_df.2$TMB_type <- factor(boxp_df.2$TMB_type, levels = c("WGS\n(3000 Mb)", "isWES\n(33 Mb)", "F1CDx assay\n(324 genes, 2.2 Mb)", "TSO500\n(523 genes, 1.94 Mb)"))

boxp_df.2 <- boxp_df.2 %>% tibble::rownames_to_column(var="outlier") %>% group_by(TMB_type) %>% mutate(is_outlier=ifelse(is_outlier(TMB), TMB, as.numeric(NA)))
boxp_df.2$outlier[which(is.na(boxp_df.2$is_outlier))] <- as.numeric(NA)
boxp_df.2$outlier[which(!is.na(boxp_df.2$is_outlier))] <- as.character(boxp_df.2$Patient_ID[!is.na(boxp_df.2$is_outlier)])

top.2 <- ggplot(boxp_df.2, aes(x=TMB_type, y=TMB, color=TMB_type)) + 
  geom_boxplot(outlier.shape = 2, outlier.size = 2.5, color = "black") +
  #ggtitle("Comparison of TMB values between sequencing methods\n(n = 93 patients)") +
  geom_jitter(shape=21, width = 0.25, size = 2) +
  geom_text_repel(data = boxp_df.2[boxp_df.2$TMB >= 15,], aes(label=outlier), na.rm = T, nudge_x = 0.25, direction = "y", force = 1.5) +
  stat_summary(fun=mean, geom="point", shape=18, size=5, color = "black") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = "dashed"), 
        panel.grid.minor.y = element_line(colour = "black", size = 0.5, linetype = "dotted"),
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        axis.line.x = element_blank(),
        plot.margin = unit(c(0.5, 0.1, 0, 0.5), "cm")) +
  scale_color_manual(values = Palette[c(1,3,5,8)], "Sequencing method") +
  scale_y_continuous(breaks = seq(10,35,5), minor_breaks = seq(7.5,35,2.5), expand=c(-0.01,0)) + 
  xlab(NULL) + 
  ylab(NULL) + 
  coord_cartesian(ylim=c(15,35.5), clip = "off")

bot.2 <- ggplot(boxp_df.2, aes(x=TMB_type, y=TMB, color=TMB_type)) + 
  geom_boxplot(notch = TRUE, outlier.shape = 2, outlier.size = 2.5, color = "black") +
  geom_jitter(shape=21, width = 0.25, size = 2) +
  geom_text_repel(data = boxp_df.2[boxp_df.2$TMB < 15,], aes(label=outlier), na.rm = T, nudge_x = 0.25, direction = "y", force = 1.5) +
  stat_summary(fun=mean, geom="point", shape=18, size=5, color = "black") +
  geom_hline(yintercept=0, size = 0.75) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = "dashed"), 
        panel.grid.minor.y = element_line(colour = "black", size = 0.5, linetype = "dotted"),
        plot.title = element_text(hjust = 0.5, size = 16),   
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14, hjust = 0.9),
        legend.position = "none",
        plot.margin = unit(c(-0.1, 0.1, 0.5, 0.5), "cm")) +
  scale_color_manual(values = Palette[c(1,3,5,8)], "Sequencing method") +
  scale_y_continuous(breaks = seq(0,15,5), minor_breaks = seq(0,15,2.5), expand=c(0.005,0)) +
  
  xlab("Sequencing method") + 
  ylab("TMB (mutations/Mb)") + 
  coord_cartesian(ylim=c(0,15))

plot_grid(top.2, bot.2, rel_heights=c(0.33, 0.66), ncol = 1, align = "v")
ggsave(file="TMB graphs/TMB_boxplot_v2.png", width=11, height=11, dpi=300)
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

###### Simple scatter plots
# WGS x isWES TMB *WITH ALL DATA POINTS*
#model.WGSxWES <- lm(NSLC_TMB$complete_WES_TMB~NSLC_TMB$complete_WGS_TMB)
ggplot(NSLC_TMB, aes(x=complete_WGS_TMB, y=complete_WES_TMB, label = Patient_ID)) + 
  geom_point() +
  stat_cor(label.x = 9, label.y = 3, color = Palette[8]) +
  stat_cor(label.x = 9, label.y = 3.8, color = Palette[8], method = "spearman", cor.coef.name = "rho") +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, label.x = 9, label.y = 2.2, color = Palette[8]) +
  geom_smooth(method = lm, se = FALSE, color = Palette[8]) +
  geom_abline(linetype="dashed", size = 1, color = Palette[3]) +
  #ggtitle("Correlation between WGS TMB and isWES TMB") +
  #geom_text(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10), ), hjust=0.5, vjust=-0.5, ) +
  #geom_label_repel(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10)), size = 3, box.padding = 0.5, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  scale_x_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  scale_y_continuous(name = "isWES TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"), plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")
  #annotate("text", x=11, y=3, label=c(paste("Pearson r = ", paste(as.character(format(cor(NSLC_TMB$complete_WGS_TMB, NSLC_TMB$complete_WES_TMB, method="pearson"), digits=3)), "\n")), "\n\n\n", 
                                      #paste("\na = ", as.character(round(model.WGSxWES$coefficients[2], digits = 3)))), color = Palette[8])

ggsave(file="TMB graphs/WGSxisWES_TMB.png", width=6, height=6, dpi=300)

# WGS x isWES TMB *with exclusion of patient 0020*
#model.WGSxWES_0020 <- lm(NSLC_TMB[!NSLC_TMB$Patient_ID=="NSLC-0020",]$complete_WES_TMB~NSLC_TMB[!NSLC_TMB$Patient_ID=="NSLC-0020",]$complete_WGS_TMB)
ggplot(NSLC_TMB[!NSLC_TMB$Patient_ID=="0020",], aes(x=complete_WGS_TMB, y=complete_WES_TMB, label = Patient_ID)) + 
  geom_point() +
  stat_cor(label.x = 6.5, label.y = 2, color = Palette[8]) +
  stat_cor(label.x = 6.5, label.y = 2.5, color = Palette[8], method = "spearman", cor.coef.name = "rho") +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, label.x = 6.5, label.y = 1.5, color = Palette[8]) +
  #ggtitle("Correlation between WGS TMB and isWES TMB\nwith the exclusion of patient 20") +
  geom_smooth(method = lm, se = FALSE, color = Palette[8]) +
  geom_abline(linetype="dashed", size = 1, color = Palette[3]) +
  #geom_text(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10), ), hjust=0.5, vjust=-0.5, ) +
  #geom_label_repel(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10)), size = 3, box.padding = 0.5, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  scale_x_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,15,1)), limits = c(NA,15), expand = c(0,0)) +
  scale_y_continuous(name = "isWES TMB (mutations/Mb)", breaks = c(seq(0,15,1)), limits = c(NA,15), expand = c(0,0)) +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"), plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")
  #annotate("text", x=7.5, y=2, label=paste("Pearson r = ", paste(as.character(format(cor(NSLC_TMB[!NSLC_TMB$Patient_ID=="NSLC-0020",]$complete_WGS_TMB, NSLC_TMB[!NSLC_TMB$Patient_ID=="NSLC-0020",]$complete_WES_TMB, method="pearson"), digits=3)), "\n"), paste("a = ", as.character(round(model.WGSxWES_0020$coefficients[2], digits = 3)))), color = "#e76f51")

ggsave(file="TMB graphs/WGSxisWES_TMB_no_NSLC-0020.png", width=6, height=6, dpi=300)

# WGS x isWES TMB *EXCLUDING ALL OUTLIERS FROM WGS AND WES*
min_outlier_WES <- min(na.omit(boxp_df$is_outlier[boxp_df$TMB_type == "isWES\n(32.8 Mb)"]))
min_outlier_WGS <- min(na.omit(boxp_df$is_outlier[boxp_df$TMB_type == "WGS\n(3000 Mb)"]))
data_no_outliers <- NSLC_TMB[NSLC_TMB$complete_WES_TMB < min_outlier_WES & NSLC_TMB$complete_WGS_TMB < min_outlier_WGS,]
#model.WGSxWES_no_outliers <- lm(data_no_outliers$complete_WES_TMB ~ data_no_outliers$complete_WGS_TMB)

ggplot(data_no_outliers, aes(x=complete_WGS_TMB, y=complete_WES_TMB, label = Patient_ID)) + 
  geom_point() +
  stat_cor(label.x = 2, label.y = 0.35, color = Palette[8]) +
  stat_cor(label.x = 2, label.y = 0.5, color = Palette[8], method = "spearman", cor.coef.name = "rho") +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, label.x = 2, label.y = 0.20, color = Palette[8]) +
  geom_smooth(method = lm, formula = y~x, se = FALSE, color = Palette[8]) +
  geom_abline(linetype="dashed", size = 1, color = Palette[3]) +
  #ggtitle("Correlation between WGS TMB and isWES TMB\nwith the exlcusion of outliers") +
  #geom_text(aes(label = substring(data_no_outliers[,"Patient_ID"], 6, 10), ), hjust=0.5, vjust=-0.5, ) +
  #geom_label_repel(aes(label = substring(data_no_outliers[,"Patient_ID"], 6, 10)), size = 3, box.padding = 0.5, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  scale_x_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,4,0.5)), limits = c(NA, 4), expand = c(0,0)) +
  scale_y_continuous(name = "isWES TMB (mutations/Mb)", breaks = c(seq(0,4,0.5)), limits = c(NA, 4), expand = c(0,0)) +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"), plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")
  #annotate("text", x=2.5, y=0.5, label=c(paste("Pearson r = ", paste(as.character(format(cor(data_no_outliers$complete_WGS_TMB, data_no_outliers$complete_WES_TMB, method="pearson"), digits=3)), "\n")), "\n\n\n", paste("\na = ", as.character(round(model.WGSxWES_no_outliers$coefficients[2], digits = 3)))), color = Palette[8])

ggsave(file="TMB graphs/WGSxisWES_TMB_no_outliers.png", width=6, height=6, dpi=300)

# WGS x isWES TMB *INCLUDING OUTLIERS FROM WGS AND WES*, like the first scatter plot, but changing outliers to a different shape.
min_outlier_WES <- min(na.omit(boxp_df$is_outlier[boxp_df$TMB_type == "isWES\n(32.8 Mb)"]))
min_outlier_WGS <- min(na.omit(boxp_df$is_outlier[boxp_df$TMB_type == "WGS\n(3000 Mb)"]))
data_no_outliers <- NSLC_TMB[NSLC_TMB$complete_WES_TMB < min_outlier_WES & NSLC_TMB$complete_WGS_TMB < min_outlier_WGS,]
data_outliers <- NSLC_TMB[NSLC_TMB$complete_WES_TMB >= min_outlier_WES | NSLC_TMB$complete_WGS_TMB >= min_outlier_WGS,]
#model.WGSxWES_outliers <- lm(NSLC_TMB$complete_WES_TMB ~ NSLC_TMB$complete_WGS_TMB)
#model.WGSxWES_no_outliers <- lm(data_no_outliers$complete_WES_TMB ~ data_no_outliers$complete_WGS_TMB)

ggplot(data_no_outliers, aes(x=complete_WGS_TMB, y=complete_WES_TMB, label = Patient_ID)) + 
  geom_point() +
  geom_point(data = data_outliers, shape=17, size = 2.2, color = "black") +
  stat_cor(label.x = 9, label.y = 1.4, color = Palette[8]) +
  stat_cor(label.x = 9, label.y = 2.1, color = Palette[8], method = "spearman", cor.coef.name = "rho") +
  stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, label.x = 9, label.y = 0.7, color = Palette[8]) +
  geom_smooth(method = lm, se = FALSE, color = Palette[8]) +
  stat_cor(data = NSLC_TMB, label.x = 9, label.y = 3.7, color = Palette[1]) +
  stat_cor(label.x = 9, label.y = 4.4, color = Palette[1], method = "spearman", cor.coef.name = "rho") +
  stat_regline_equation(data = NSLC_TMB, aes(label = ..eq.label..), formula = y ~ x, label.x = 9, label.y = 3, color = Palette[1]) +
  geom_smooth(data = NSLC_TMB, method = lm, se = FALSE, color = Palette[1]) +
  geom_abline(linetype="dashed", size = 1, color = Palette[3]) +
  #ggtitle("Correlation between WGS TMB and isWES TMB") +
  #geom_text(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10), ), hjust=0.5, vjust=-0.5, ) +
  #geom_label_repel(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10)), size = 3, box.padding = 0.5, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  scale_x_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  scale_y_continuous(name = "isWES TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"), plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")
  #annotate("text", x=11, y=3, label=c(paste("Pearson r = ", paste(as.character(format(cor(NSLC_TMB$complete_WGS_TMB, NSLC_TMB$complete_WES_TMB, method="pearson"), digits=3)))),
                                      #paste("\n\na = ", as.character(round(model.WGSxWES_outliers$coefficients[2], digits = 3))),
                                      #paste("\n\n\n\nPearson r = ", paste(as.character(format(cor(data_no_outliers$complete_WGS_TMB, data_no_outliers$complete_WES_TMB, method="pearson"), digits=3)))), 
                                      #paste("\n\n\n\n\n\na = ", as.character(round(model.WGSxWES_no_outliers$coefficients[2], digits = 3)))),
           #color = c(Palette[6], Palette[6], Palette[8], Palette[8]))

ggsave(file="TMB graphs/WGSxisWES_TMB_outliers.png", width=6, height=6, dpi=300)

# WGS x F1CDx NGS TMB
ggplot(NSLC_TMB, aes(x=complete_WGS_TMB, y=F1CDx_TMB, label = Patient_ID)) + 
  geom_point() +
  #ggtitle("Correlation between WGS TMB and F1CDx panel TMB") +
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE) +
  geom_abline(linetype="dashed", size = 1, color = Palette[3]) +
  #geom_text(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10), ), hjust=0.5, vjust=-0.5, ) +
  #geom_label_repel(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10)), size = 3, box.padding = 0.5, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  scale_x_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  scale_y_continuous(name = "F1CDx TMB (mutations/Mb)", breaks = c(seq(0,35,1)), limits = c(NA, 35), expand = c(0,0)) +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/WGSxF1CDx_TMB.png", width=5, height=5, dpi=300)

# isWES x F1CDx NGS TMB
ggplot(NSLC_TMB, aes(x=complete_WES_TMB, y=F1CDx_TMB, label = Patient_ID)) + 
  geom_point() +
  #ggtitle("Correlation between isWES TMB and F1CDx TMB") +
  geom_abline(linetype="dashed", size = 1, color = Palette[3]) +
  #geom_smooth(method = lm, se = FALSE, color = "red") +
  #geom_text(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10), ), hjust=0.5, vjust=-0.5, ) +
  #geom_label_repel(aes(label = substring(NSLC_TMB[,"Patient_ID"], 6, 10)), size = 3, box.padding = 0.5, point.padding = 1, max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  scale_x_continuous(name = "isWES TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0,0)) +
  scale_y_continuous(name = "F1CDx TMB (mutations/Mb)", breaks = c(seq(0,35,1)), limits = c(NA, 35), expand = c(0,0)) +
  theme(plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/isWESxF1CDx_TMB.png", width=5, height=5, dpi=300)


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

###### Bland-Altman analysis of pairs of TMB values (inspired by Heeke et al. [2020], figure 2, https://doi.org/10.1016/j.jtho.2020.05.013)
# WGS x isWES
bdr_WXS.stats <- blandr.statistics(NSLC_TMB$complete_WGS_TMB, NSLC_TMB$complete_WES_TMB)
bdr_WXS.diff_mean <- round(mean(bdr_WXS.stats$differences), digits = 2)
blandr.draw(NSLC_TMB$complete_WGS_TMB, NSLC_TMB$complete_WES_TMB, method1name = "WGS", method2name = "isWES", ciDisplay = F, ciShading = F, plotTitle = "Bland-Altman plot for comparison of WGS and isWES") +
  geom_errorbar(aes(x=7.5, ymin=bdr_WXS.stats$biasLowerCI, ymax=bdr_WXS.stats$biasUpperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=7.5, ymin=bdr_WXS.stats$upperLOA_lowerCI, ymax=bdr_WXS.stats$upperLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=7.5, ymin=bdr_WXS.stats$lowerLOA_lowerCI, ymax=bdr_WXS.stats$lowerLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_hline(aes(yintercept=bdr_WXS.stats$bias), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_WXS.stats$upperLOA), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_WXS.stats$lowerLOA), linetype = "dashed", color = Palette[5]) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = c("dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "blank")),
        axis.text.y = element_text(color = c("black", "black", "black", "black", "black", "black", "black", "black", Palette[5]))) +
  scale_x_continuous(limits = c(0,15), expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq(-2,3,0.5), bdr_WXS.diff_mean)) +
  xlab("TMB means (mutations/Mb)") +
  ylab("TMB differences (mutations/Mb)") +
  coord_cartesian(clip = "off")
  
ggsave(file="TMB graphs/Bland-Altman plots/BA_WGSxisWES_TMB.png", width=10, height=10, dpi=300)

# WGS x Illu500 NGS
bdr_WGS_Illu.stats <- blandr.statistics(NSLC_TMB$complete_WGS_TMB, NSLC_TMB$Illu500_TMB)
bdr_WGS_Illu.diff_mean <- round(mean(bdr_WGS_Illu.stats$differences), digits = 2)
blandr.draw(NSLC_TMB$complete_WGS_TMB, NSLC_TMB$Illu500_TMB, method1name = "WGS", method2name = "Illu500", ciDisplay = F, ciShading = F, plotTitle = "Bland-Altman plot for comparison of WGS and Illumina TruSight Oncology 500 panel") +
  geom_errorbar(aes(x=10, ymin=bdr_WGS_Illu.stats$biasLowerCI, ymax=bdr_WGS_Illu.stats$biasUpperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=10, ymin=bdr_WGS_Illu.stats$upperLOA_lowerCI, ymax=bdr_WGS_Illu.stats$upperLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=10, ymin=bdr_WGS_Illu.stats$lowerLOA_lowerCI, ymax=bdr_WGS_Illu.stats$lowerLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_hline(aes(yintercept=bdr_WGS_Illu.stats$bias), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_WGS_Illu.stats$upperLOA), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_WGS_Illu.stats$lowerLOA), linetype = "dashed", color = Palette[5]) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = c("dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "blank")),
        axis.text.y = element_text(color = c("black", "black", "black", "black", "black", "black", Palette[5]))) +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq(-20,10,5), bdr_WGS_Illu.diff_mean), limits = c(-15,10)) +
  xlab("TMB means (mutations/Mb)") +
  ylab("TMB differences (mutations/Mb)") +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/Bland-Altman plots/BA_WGSxIllu500_TMB.png", width=10, height=10, dpi=300)

# isWES x Illu500 NGS
bdr_WES_Illu.stats <- blandr.statistics(NSLC_TMB$complete_WES_TMB, NSLC_TMB$Illu500_TMB)
bdr_WES_Illu.diff_mean <- round(mean(bdr_WES_Illu.stats$differences), digits = 2)
blandr.draw(NSLC_TMB$complete_WES_TMB, NSLC_TMB$Illu500_TMB, method1name = "isWES", method2name = "Illu500", ciDisplay = F, ciShading = F, plotTitle = "Bland-Altman plot for comparison of isWES and Illumina TruSight Oncology 500 panel") +
  geom_errorbar(aes(x=10, ymin=bdr_WES_Illu.stats$biasLowerCI, ymax=bdr_WES_Illu.stats$biasUpperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=10, ymin=bdr_WES_Illu.stats$upperLOA_lowerCI, ymax=bdr_WES_Illu.stats$upperLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=10, ymin=bdr_WES_Illu.stats$lowerLOA_lowerCI, ymax=bdr_WES_Illu.stats$lowerLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_hline(aes(yintercept=bdr_WES_Illu.stats$bias), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_WES_Illu.stats$upperLOA), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_WES_Illu.stats$lowerLOA), linetype = "dashed", color = Palette[5]) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = c("dotted", "dotted", "dotted", "blank")),
        axis.text.y = element_text(color = c("black", "black", "black", "black", "black", Palette[5]))) +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq(-20,10,5), bdr_WES_Illu.diff_mean), limits = c(-15,5)) +
  xlab("TMB means (mutations/Mb)") +
  ylab("TMB differences (mutations/Mb)") +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/Bland-Altman plots/BA_isWESxIllu500_TMB.png", width=10, height=10, dpi=300)

# Illu500 NGS x F1CDx NGS
bdr_Illu_F1CDx.stats <- blandr.statistics(NSLC_TMB$Illu500_TMB, NSLC_TMB$F1CDx_TMB)
bdr_Illu_F1CDx.diff_mean <- round(mean(bdr_Illu_F1CDx.stats$differences), digits = 2)
blandr.draw(NSLC_TMB$Illu500_TMB, NSLC_TMB$F1CDx_TMB, method1name = "Illu500", method2name = "F1CDx", ciDisplay = F, ciShading = F, plotTitle = "Bland-Altman plot for comparison of Illumina TruSight Oncology 500 panel\nand Foundation Medicine FoundationOneCDx panel") +
  geom_errorbar(aes(x=10, ymin=bdr_Illu_F1CDx.stats$biasLowerCI, ymax=bdr_Illu_F1CDx.stats$biasUpperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=10, ymin=bdr_Illu_F1CDx.stats$upperLOA_lowerCI, ymax=bdr_Illu_F1CDx.stats$upperLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_errorbar(aes(x=10, ymin=bdr_Illu_F1CDx.stats$lowerLOA_lowerCI, ymax=bdr_Illu_F1CDx.stats$lowerLOA_upperCI), width = 0.1, color = Palette[3]) +
  geom_hline(aes(yintercept=bdr_Illu_F1CDx.stats$bias), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_Illu_F1CDx.stats$upperLOA), linetype = "dashed", color = Palette[5]) +
  geom_hline(aes(yintercept=bdr_Illu_F1CDx.stats$lowerLOA), linetype = "dashed", color = Palette[5]) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(colour = "black", size = 0.5, linetype = c("dotted", "dotted", "dotted", "dotted", "blank")),
        axis.text.y = element_text(color = c("black", "black", "black", "black", Palette[5]))) +
  scale_x_continuous(limits = c(0,20), expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq(-10,5,5), bdr_Illu_F1CDx.diff_mean), limits = c(-10,5)) +
  xlab("TMB means (mutations/Mb)") +
  ylab("TMB differences (mutations/Mb)") +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/Bland-Altman plots/BA_Illu500xF1CDx_TMB.png", width=10, height=10, dpi=300)

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

###### Connected scatter plots
colourCount <- length(unique(boxp_df$Patient_ID))
Palette2 <- colorRampPalette(Palette)(colourCount)

# Connected plot for only WGS and isWES
csWXS_df <- rbind(boxp_df[boxp_df$TMB_type == "isWES\n(32.8 Mb)",], boxp_df[boxp_df$TMB_type == "WGS\n(3000 Mb)",])[c(2,3,4)]
csWXS_df["difference"] <- rep(ifelse(boxp_df[boxp_df$TMB_type == "WGS\n(3000 Mb)","TMB"] > boxp_df[boxp_df$TMB_type == "isWES\n(32.8 Mb)","TMB"], "WGS > isWES", "WGS < isWES"), 2)
ggplot(csWXS_df, aes(x=TMB_type, y=TMB, group = Patient_ID, color = difference)) +
  #ggtitle("Comparison of TMB values between WGS and isWES") +
  geom_point(data = csWXS_df[csWXS_df$difference == "WGS > isWES",]) +
  geom_line(data = csWXS_df[csWXS_df$difference == "WGS > isWES",], size = 0.3) +
  geom_point(data = csWXS_df[csWXS_df$difference == "WGS < isWES",]) +
  geom_line(data = csWXS_df[csWXS_df$difference == "WGS < isWES",], size = 0.7) +
  scale_color_manual(values = Palette[c(8,3)]) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.95, 0.5), plot.margin=unit(c(5.5,40,5.5,5.5),"pt")) +
  scale_y_log10(name = "TMB (mutations/Mb)", limits = c(NA,20)) +
  scale_x_discrete(name = "Sequencing method") +
  coord_cartesian(clip = "off") + 
  labs(color="TMB comparison")

ggsave(file="TMB graphs/TMB_connected_scatter_WGS_WES.png", width=6, height=6, dpi=300)

nrow(csWXS_df[csWXS_df$difference=="WGS > isWES",])/2 # nb of patient with higher WGS TMB than isWES TMB
nrow(csWXS_df[csWXS_df$difference=="WGS < isWES",])/2 # nb of patient with lower WGS TMB than isWES TMB

# Connected plot for all sequencing methods
ggplot(boxp_df, aes(x=TMB_type, y=TMB, group = as.character(Patient_ID), color = as.character(Patient_ID))) +
  ggtitle("Comparison of TMB values between WGS,\nisWES and NGS genes panels (n = 93 patients)") +
  geom_point() + 
  geom_line(size = 0.2) +
  geom_text_repel(aes(label=outlier), na.rm = T, direction = "y") +
  scale_color_manual(values = Palette2) +
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(name = "TMB (mutations/Mb)") +
  scale_x_discrete(name = "Sequencing method") +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/TMB_connected_scatter_all.png", width=10, height=10, dpi=300)


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### TMB plot where patients are grouped accroding to their histology
# Dataset for the TMB plot. Complete records of patients clinical characteristics and TMB values. *Obtained from the script "13-TMB_cutoffs.R"*
histology_df.TMB <- NSLC_TMB[c("Patient_ID", "histology", "complete_WGS_TMB")]
histology_df.TMB$Patient_ID<-gsub("NSLC-","",as.character(histology_df.TMB$Patient_ID))
histology_df.TMB$histology = factor(histology_df.TMB$histology, levels=c('3','2','4','6','7'))
# Corresponding levels : 3 = adenocarcinoma, 2 = squamous cell carcinoma, 4 = carcinoid tumor, 6 = adenosquamous carcinoma, 7 = sarcomatoid carcinoma
levels(histology_df.TMB$histology) <- c('adenocarcinoma', 'squamous cell\ncarcinoma', 'carcinoid\ntumor', 'adenosquamous\ncarcinoma', 'sarcomatoid\ncarcinoma')

# Histology plot with TMB (y axis)
png("TMB graphs/TMB_histology.png", width = 12, height = 6, units = "in", res = 1200)
{histology_plot <- ggplot(histology_df.TMB, aes(x=reorder(Patient_ID, -complete_WGS_TMB), y=complete_WGS_TMB)) +
    #ggtitle("Distribution of histological types TMB values") +
    geom_col(width=0.9, fill = Palette[5]) +
    theme_classic() +
    facet_grid(~histology, scales = "free_x", space = "free_x", switch = "x") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6, colour = "#264653"), axis.text.y = element_text(colour = "#264653"), axis.ticks.x=element_blank(),
          panel.grid.major.y = element_line(colour="black", size=0.5), legend.text=element_text(size=11), plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(), panel.spacing = unit(0.5, "lines"), strip.placement = "outside", strip.text.x = element_text(angle = 90, colour = "#264653", size = 10), 
          strip.switch.pad.grid = unit(-1, "cm"), plot.margin=unit(c(5.5,5.5,35,5.5),"pt"), panel.grid.minor = element_line("#2a9d8f", linetype = "dotted")) +
    scale_y_continuous(name="WGS TMB (mutations/Mb)",  breaks = c(seq(0,20,5)), minor_breaks = c(seq(0,20,1)), limits = c(NA, 20), expand = c(0, 0)) +
    scale_x_discrete(name="Patients")

  gp <- ggplotGrob(histology_plot)
  
  t_title <- gp$layout[gp$layout[['name']] == 'xlab-b' ,][['b']]
  t_strip <- gp$layout[grepl('strip', gp$layout[['name']]),][['b']]
  
  gp$layout[gp$layout[['name']] == 'xlab-b' ,][['t']] <- unique(t_strip)
  gp$layout[gp$layout[['name']] == 'xlab-b' ,][['b']] <- unique(t_strip)
  gp$layout[grepl('strip', gp$layout[['name']]),][['t']] <- t_title
  gp$layout[grepl('strip', gp$layout[['name']]),][['b']] <- t_title
  
  for(i in which(grepl("strip-b", gp$layout$name))){
    gp$grobs[[i]]$layout$clip <- "off"
  }
  
  for (k in grep("strip-b",gp$layout$name)) {
    gp$grobs[[k]]$grobs[[1]]$children[[1]] <- grid.lines(x=unit(c(0,1),"npc"), y=unit(c(1,1),"npc"), 
                                                         gp=gpar(col="#2a9d8f", lwd=2))
  }
  
  grid.newpage()
  grid.draw(gp)
}
dev.off()


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### TMB plot showing comparison of WGS, WES and NGS panels TMBs
TMB.comp <- NSLC_TMB[order(-NSLC_TMB$complete_WGS_TMB),]
TMB.comp <- reshape2::melt(TMB.comp[, c("Patient_ID", "complete_WGS_TMB", "complete_WES_TMB", "F1CDx_TMB", "Illu500_TMB", "MSKI468_TMB", "OncoMineTMB_TMB", "QIAseq_TMB")], id.vars="Patient_ID", variable.name = "TMB type", value.name = "TMB")
TMB.comp$Patient_ID<-gsub("NSLC-","",as.character(TMB.comp$Patient_ID))
#TMB.comp$`TMB type` <- factor(TMB.comp$`TMB type`, levels = c("WGS", "isWES", "F1CDx assay", "TSO500", "MSK-IMPACT", "NEOPlus RUO assay", "Oncomine TML assay", "QIAseq TMB panel"))
TMB.comp$`TMB type` <- as.character(TMB.comp$`TMB type`)
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "complete_WGS_TMB"] <- "WGS"
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "complete_WES_TMB"] <- "isWES"
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "F1CDx_TMB"] <- "F1CDx"
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "Illu500_TMB"] <- "TSO500"
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "MSKI468_TMB"] <- "MSK-IMPACT"
#TMB.comp$`TMB type`[TMB.comp$`TMB type` == "NeoPlus_TMB"] <- "NEOPlus RUO"
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "OncoMineTMB_TMB"] <- "Oncomine TML"
TMB.comp$`TMB type`[TMB.comp$`TMB type` == "QIAseq_TMB"] <- "QIAseq"

ggplot() +
  geom_col(data = TMB.comp[TMB.comp$`TMB type` == "WGS",], aes(x=fct_inorder(Patient_ID), y=TMB), fill = Palette[5]) +
  geom_point(data = TMB.comp[TMB.comp$`TMB type` != "WGS",], 
             aes(x=Patient_ID, y=TMB, fill = fct_inorder(`TMB type`),
                 shape = fct_inorder(`TMB type`), size = fct_inorder(`TMB type`), color = fct_inorder(`TMB type`))) +
  geom_hline(yintercept = 1.70, color = Palette[8], size = 1, linetype = "dashed") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, colour = "#264653"), axis.text.y = element_text(colour = "#264653", size = 12),
        axis.ticks.x=element_blank(), panel.grid.major.y = element_line(colour="black", size=0.5), legend.text=element_text(size=11), 
        panel.spacing = unit(0.5, "lines"), strip.placement = "outside", strip.text.x = element_text(angle = 90, colour = "#264653", size = 10), 
        plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"), panel.grid.minor = element_line("#2a9d8f", linetype = "dotted"), axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18)) +
  scale_y_continuous(name="TMB (mutations/Mb)",  breaks = c(seq(0,35,5)), minor_breaks = c(seq(0,35,1)), limits = c(NA, 35), expand = c(0, 0)) +
  scale_x_discrete(name="Patients") +
  scale_fill_manual(values=c(NA, "black", rep(NA, times = 4)), name = "TMB types\n") + # In order of legend (top to bottom). isWES "-" symbol has no fill possible.
  scale_shape_manual(values=c(45, 21, 2, 5, 0, 1), name = "TMB types\n") + 
  scale_size_manual(values=c(11.5, 2, 2, 2, 2, 2), name = "TMB types\n") +
  scale_color_manual(values=c(Palette[3], rep("black", times = 5)), name = "TMB types\n") + 
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/TMB_comparison.png", width=16, height=8, dpi=300)

# Same graph, but with a zoom on TMB ranging from 0 to 5 mut/Mb.
comp.zoom.top <- ggplot() +
  geom_col(data = TMB.comp[TMB.comp$`TMB type` == "WGS",], aes(x=fct_inorder(Patient_ID), y=TMB), fill = Palette[5], width = 0.8) +
  geom_point(data = TMB.comp[TMB.comp$`TMB type` != "WGS" & TMB.comp$TMB >=5,], 
             aes(x=Patient_ID, y=TMB, fill = fct_inorder(`TMB type`),
                 shape = fct_inorder(`TMB type`), size = fct_inorder(`TMB type`), color = fct_inorder(`TMB type`))) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "#264653", size = 12),
        axis.ticks.x = element_blank(), 
        panel.grid.major.y = element_line(colour="black", size=0.5), 
        legend.position = "none",
        panel.spacing = unit(0.5, "lines"), 
        strip.placement = "outside", 
        strip.text.x = element_text(angle = 90, colour = "#264653", size = 10), 
        plot.margin=unit(c(5,5,2,5),"pt"), 
        panel.grid.minor = element_line("#2a9d8f", linetype = "dotted"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_discrete(expand = c(0.015,0)) +
  scale_y_continuous(name="TMB (mutations/Mb)", expand = c(0, 0), breaks = c(seq(5,35,5)), minor_breaks = c(seq(5,35,1))) +
  scale_fill_manual(values=c(NA, "black", rep(NA, times = 4)), name = "TMB types\n") + # In order of legend (top to bottom). isWES "-" symbol has no fill possible.
  scale_shape_manual(values=c(45, 21, 2, 5, 0, 1), name = "TMB types\n") + 
  scale_size_manual(values=c(11.5, 2, 2, 2, 2, 2), name = "TMB types\n") +
  scale_color_manual(values=c(Palette[3], rep("black", times = 5)), name = "TMB types\n") + 
  coord_cartesian(ylim = c(5,35), clip = "off")

comp.zoom.bot <- ggplot() +
  geom_col(data = TMB.comp[TMB.comp$`TMB type` == "WGS",], aes(x=fct_inorder(Patient_ID), y=TMB), fill = Palette[5], width = 0.8) +
  geom_point(data = TMB.comp[TMB.comp$`TMB type` != "WGS" & TMB.comp$TMB<5,], 
             aes(x=Patient_ID, y=TMB, fill = fct_inorder(`TMB type`),
                 shape = fct_inorder(`TMB type`), size = fct_inorder(`TMB type`), color = fct_inorder(`TMB type`))) +
  geom_hline(yintercept = 1.70, color = Palette[8], size = 1, linetype = "dashed") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, colour = "#264653"), 
        axis.text.y = element_text(colour = "#264653", size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(colour="black", size=0.5), 
        legend.position = "none", 
        panel.spacing = unit(0.5, "lines"), 
        strip.placement = "outside", 
        strip.text.x = element_text(angle = 90, colour = "#264653", size = 10), 
        plot.margin=unit(c(-1,5,5,5),"pt"), 
        panel.grid.minor = element_line("#2a9d8f", linetype = "dotted"), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18, hjust = 1.5),
        axis.line.x = element_blank()) +
  scale_y_continuous(name="TMB (mutations/Mb)", expand = c(0,0.05), breaks = c(seq(0,5,5)), minor_breaks = c(seq(0,5,1))) +
  scale_x_discrete(name="Patients", expand = c(0.015,0)) +
  scale_fill_manual(values=c(NA, "black", rep(NA, times = 4)), name = "TMB types\n") + # In order of legend (top to bottom). isWES "-" symbol has no fill possible.
  scale_shape_manual(values=c(45, 21, 2, 5, 0, 1), name = "TMB types\n") + 
  scale_size_manual(values=c(11.5, 2, 2, 2, 2, 2), name = "TMB types\n") +
  scale_color_manual(values=c(Palette[3], rep("black", times = 5)), name = "TMB types\n") + 
  coord_cartesian(ylim = c(0,4.9), clip = "on")

plot_grid(comp.zoom.top, comp.zoom.bot, rel_heights=c(0.40, 0.60), ncol = 1, align = "v")
ggsave(file="TMB graphs/TMB_comparison_zoom.png", width=16, height=8, dpi=300)


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Clinicopathological plots

## Agreement of WGS and histology plot
# Dataframe
cp.WGSxhisto.df <- NSLC_TMB[,c("histology", "complete_WGS_TMB")]
cp.WGSxhisto.df$histology <- factor(cp.WGSxhisto.df$histology, levels = c(2,3,6,7), labels = c("Squamous cell\ncarcinoma\nn=2", 
                                                                                               "Adenocarcinoma\nn=83", 
                                                                                               "Adenosquamous\ncarcinoma\nn=2", 
                                                                                               "Sarcomatoid\ncarcinoma\nn=6"))

cp.WGSxhisto.df["outlier"] <- NA
cp.WGSxhisto.outliers <- lapply(unique(cp.WGSxhisto.df$histology), 
                   FUN = function(x) is_outlier(cp.WGSxhisto.df$complete_WGS_TMB[cp.WGSxhisto.df$histology == x]))
cp.WGSxhisto.df$outlier[cp.WGSxhisto.df$histology == unique(cp.WGSxhisto.df$histology)[1]] <- cp.WGSxhisto.outliers[[1]]
cp.WGSxhisto.df$outlier[cp.WGSxhisto.df$histology == unique(cp.WGSxhisto.df$histology)[2]] <- cp.WGSxhisto.outliers[[2]]
cp.WGSxhisto.df$outlier[cp.WGSxhisto.df$histology == unique(cp.WGSxhisto.df$histology)[3]] <- cp.WGSxhisto.outliers[[3]]
cp.WGSxhisto.df$outlier[cp.WGSxhisto.df$histology == unique(cp.WGSxhisto.df$histology)[4]] <- cp.WGSxhisto.outliers[[4]]

# P-value
cp.WGSxhisto.p <- round(summary(aov(complete_WGS_TMB~histology, data=NSLC_TMB))[[1]]$Pr[1], digits = 3)

# Plot
ggplot(cp.WGSxhisto.df, aes(x=reorder(histology, -complete_WGS_TMB), y=complete_WGS_TMB)) +
  geom_boxplot(outlier.shape = 2, outlier.size = 2.5, color = "black", width = 0.5) +
  geom_jitter(data = cp.WGSxhisto.df[!cp.WGSxhisto.df$outlier,], position=position_jitter(0.15), shape=21, cex=1.3, fill = Palette[5]) +
  annotate(geom = "text", label = paste0("ANOVA p-value = ",cp.WGSxhisto.p), x=4, y=16, size = 4) +
  theme_classic() +
  scale_y_continuous(name = "WGS TMB (mutations/Mb)", limits = c(0,20), breaks = c(seq(0,20,1)), expand = c(0,0)) +
  scale_x_discrete(name = "Histological types") +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"),
        plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/cp_WGSxhistology.png", width=6, height=6, dpi=300)
  

## Agreement of WGS and pathological stage plot
# Dataframe
cp.WGSxpath.df <- NSLC_TMB[c("path.stage","complete_WGS_TMB")]
cp.WGSxpath.df["path.stage"] <- ifelse(NSLC_TMB$path.stage %in% c("1A1","1A2","1A3","1B"), "Stage I", ifelse(NSLC_TMB$path.stage %in% c("2A", "2B"), "Stage II", "Stage III & IV"))

cp.WGSxpath.group <- aggregate(cp.WGSxpath.df[,"path.stage"], list(cp.WGSxpath.df$path.stage), length)
colnames(cp.WGSxpath.group) <- c("path.stage", "N")
cp.WGSxpath.group["max"] <- aggregate(cp.WGSxpath.df[,c("complete_WGS_TMB","path.stage")], list(cp.WGSxpath.df$path.stage), max)[,2]

# P-value
cp.WGSxpath.p <- round(summary(aov(complete_WGS_TMB~path.stage, data=NSLC_TMB))[[1]]$Pr[1], digits = 3)

# Plot. Reference : http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
ggplot(cp.WGSxpath.df, aes(x=path.stage, y=complete_WGS_TMB, group = path.stage)) +
  geom_violin(fill=Palette[5]) +
  geom_boxplot(width=0.05)+
  geom_text(data=cp.WGSxpath.group, aes(label=paste0("n = ",N), y=max, fontface = 2), nudge_y = 0.25) +
  annotate(geom = "text", label = paste0("ANOVA p-value < ",cp.WGSxpath.p), x=1, y=16, size = 4) +
  theme_classic() +
  scale_y_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(0, 20), expand = c(0.01,0)) +
  scale_x_discrete(name = "Pathological stages") +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"),
        plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/cp_WGSxpath_stage.png", width=6, height=6, dpi=300)


## Agreement of WGS and age+sex plot
# Dataframe
cp.WGSxage_sex.df <- NSLC_TMB[,c("Patient_ID", "sex", "age", "complete_WGS_TMB")]
cp.WGSxage_sex.df$sex <- ifelse(cp.WGSxage_sex.df$sex==2, "Female", "Male")
cp.WGSxage_sex.df$age_group <- ifelse(cp.WGSxage_sex.df$age<60, "<60", ifelse(cp.WGSxage_sex.df$age>=60 & cp.WGSxage_sex.df$age<70, "60-70", ">=70"))
cp.WGSxage_sex.df$age_group <- factor(cp.WGSxage_sex.df$age_group, levels = c("<60", "60-70", ">=70"))

# P-values
cp.WGSxage_sex.age_p <- round(summary(aov(complete_WGS_TMB~age_group, data=cp.WGSxage_sex.df))[[1]]$Pr[1], digits = 3)

var(NSLC_TMB[NSLC_TMB$sex==1,]$complete_WGS_TMB)/var(NSLC_TMB[NSLC_TMB$sex==2,]$complete_WGS_TMB) # Ratio of variances is >4, therefore the sex variances are not equal
cp.WGSxage_sex.sex_p <- round(t.test(complete_WGS_TMB~sex, data = NSLC_TMB, var.equal = FALSE)$p.value, digits = 3)

# Dataframe for p-values position in facet_grid(). Reference: https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
ann_text_age_sex <- data.frame(Patient_ID = "NSLC-0007", complete_WGS_TMB = 13.5, sex = "Female",
                       age_group = factor(">=70", levels = c("<60", "60-70", ">=70")))

# Plot
png("TMB graphs/cp_WGSxage_sex.png", width = 6, height = 6, units = "in", res = 1200)
{cp.WGSxage_sex.plot <- ggplot(cp.WGSxage_sex.df, aes(x=reorder(Patient_ID, -complete_WGS_TMB), y=complete_WGS_TMB, fill = sex)) +
    geom_col(width = 0.7) +
    facet_grid(~age_group, scales = "free", space = "free", shrink = FALSE, as.table = FALSE, switch = "x", 
               labeller = as_labeller(c("<60" = "<60\nn=22",
                                        "60-70" = "60-70\nn=37",
                                        ">=70" = "\u226570\nn=34"))) +
    #geom_text(data=ann_text_age_sex, label = paste0("Sex t-test p-value = ", cp.WGSxage_sex.sex_p), hjust = 0, nudge_y = -0.5, nudge_x = -25, size = 4) +
    #geom_text(data=ann_text_age_sex, label = paste0("\n\nAge ANOVA p-value = ", cp.WGSxage_sex.age_p), hjust = 0, nudge_y = -0.5, nudge_x = -25, size = 4) +
    theme_classic() +
    scale_y_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(0, 20), expand = c(0,0)) +
    scale_x_discrete(name="Age") +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(colour = "#264653"),
          plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5), axis.ticks.x = element_blank(),
          legend.position=c(.9,.75), axis.line.x.bottom = element_blank(), strip.text.x = element_text(colour = "#264653", size = 10),
          strip.switch.pad.grid = unit(0, "cm"), strip.placement = "outside") +
    scale_fill_manual(values = c(Palette[c(5,3)]), name = "Sex") +
    coord_cartesian(clip = "off")

  gp <- ggplotGrob(cp.WGSxage_sex.plot)

  for(i in which(grepl("strip-b", gp$layout$name))){
    gp$grobs[[i]]$layout$clip <- "off"
  }
  
  for (k in grep("strip-b",gp$layout$name)) {
    gp$grobs[[k]]$grobs[[1]]$children[[1]] <- grid.lines(x=unit(c(0,1),"npc"), y=unit(c(1,1),"npc"), 
                                                         gp=gpar(col=Palette[8], lwd=2))
  }
  
  grid.newpage()
  grid.draw(gp)
  }
dev.off()


## Agreement of WGS and tumor size plot
# Dataframe
cp.WGSxsize.df <- NSLC_TMB[,c("tumor.size (cm)", "complete_WGS_TMB", "path.stage")]
cp.WGSxsize.df$`tumor.size (cm)` <- round(as.numeric(cp.WGSxsize.df$`tumor.size (cm)`), digits = 2)
cp.WGSxsize.df <- cp.WGSxsize.df %>% mutate(size_group = case_when(`tumor.size (cm)` <=3 ~ "<=3",
                                                 `tumor.size (cm)` > 3 & `tumor.size (cm)` <=5 ~ ">3 - <=5",
                                                 `tumor.size (cm)` > 5 & `tumor.size (cm)` <=7 ~ ">5 - <=7",
                                                 `tumor.size (cm)` > 7 ~ ">7"))
cp.WGSxsize.df["outlier"] <- NA
cp.WGSxsize.outliers <- lapply(unique(cp.WGSxsize.df$size_group), 
                   FUN = function(x) is_outlier(cp.WGSxsize.df$complete_WGS_TMB[cp.WGSxsize.df$size_group == x]))
cp.WGSxsize.df$outlier[cp.WGSxsize.df$size_group == unique(cp.WGSxsize.df$size_group)[1]] <- cp.WGSxsize.outliers[[1]]
cp.WGSxsize.df$outlier[cp.WGSxsize.df$size_group == unique(cp.WGSxsize.df$size_group)[2]] <- cp.WGSxsize.outliers[[2]]
cp.WGSxsize.df$outlier[cp.WGSxsize.df$size_group == unique(cp.WGSxsize.df$size_group)[3]] <- cp.WGSxsize.outliers[[3]]
cp.WGSxsize.df$outlier[cp.WGSxsize.df$size_group == unique(cp.WGSxsize.df$size_group)[4]] <- cp.WGSxsize.outliers[[4]]

cp.WGSxsize.group <- aggregate(cp.WGSxsize.df[,"size_group"], list(cp.WGSxsize.df$size_group), length)
cp.WGSxsize.group["mean_WGS_TMB"] <- aggregate(cp.WGSxsize.df[,"complete_WGS_TMB"], list(cp.WGSxsize.df$size_group), mean)[,2]
colnames(cp.WGSxsize.group) <- c("size_group", "n", "mean_WGS_TMB")

# P-value
cp.WGSxsize.p <- round(summary(aov(complete_WGS_TMB~size_group, data=cp.WGSxsize.df))[[1]]$Pr[1], digits = 3)


# Plot
ggplot(cp.WGSxsize.df, aes(x=size_group, y=complete_WGS_TMB)) +
  geom_boxplot(outlier.shape = 2, outlier.size = 2.5, color = "black", width = 0.5) +
  geom_jitter(data = cp.WGSxsize.df[!cp.WGSxsize.df$outlier,], position=position_jitter(0.15), shape=21, cex=1.3, fill = Palette[5]) +
  annotate("text", label = paste0("ANOVA p-value < ", cp.WGSxsize.p), hjust = 0, x=0.5, y=16, size = 4) +
  theme_classic() +
  scale_y_continuous(name = "WGS TMB (mutations/Mb)", breaks = c(seq(0,20,1)), limits = c(0, 20), expand = c(0,0)) +
  scale_x_discrete(name="Tumor size (cm)", labels = c("<=3" = "\u22643\nn=58",
                                                      ">3 - <=5" = ">3 - \u22645\nn=19",
                                                      ">5 - <=7" = ">5 - \u22647\nn=12",
                                                      ">7" = ">7\nn=4")) +
  theme(axis.text.x = element_text(colour = "#264653"), axis.text.y = element_text(colour = "#264653"),
        plot.margin=unit(c(5.5,10,5.5,5.5),"pt"), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(clip = "off")

ggsave(file="TMB graphs/cp_WGSxsize.png", width=6, height=6, dpi=300)


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# Other infos

# Mutations counts info (WGS and WES)
mean(NSLC_TMB$WGS_mutation_count)/mean(NSLC_TMB$WES_mutation_count)
median(NSLC_TMB$WGS_mutation_count)/median(NSLC_TMB$WES_mutation_count)

# Panels 0 mutations/Mb patients count
nrow(NSLC_TMB[NSLC_TMB$F1CDx_TMB == 0,])
nrow(NSLC_TMB[NSLC_TMB$Illu500_TMB == 0,])



