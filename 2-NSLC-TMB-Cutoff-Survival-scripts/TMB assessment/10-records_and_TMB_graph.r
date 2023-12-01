# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

library(ggplot2)
library(data.table)
library(grid)
library(tidyr)
library(gridExtra)

# RStudio execution
#output_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#message(paste("Output directory :", output_dir))

# Command line execution (Rscript)
#setwd("/mnt/raid5_hdd/ruelou01")
#message(paste("Output directory: ", "/mnt/raid5_hdd/ruelou01", sep=""))

patients.info <- readxl::read_excel("NCI_NeverSmoker_n131_20210812_share_updated-124_125_fix.xlsx", range = "A1:W130")
patients.id <- readxl::read_excel("Table_ID-124_125_fix.xlsx", range = "A1:W125")[,c(1,5)]

records <- read.csv2("coding_non-coding_records.csv")
records["Sample_ID"] <- merge(records, patients.id, by.x = "Patient_ID", by.y = "Subject")[,"Tumor_Sample_Name"]
records <- merge(records, patients.info[,c("ID","histology")], by.x = "Sample_ID", by.y = "ID", sort = FALSE, all = TRUE)[,c(2,1,3:7)]
records <- records[!is.na(records$Patient_ID),]
records <- records[order(records$Patient_ID),]
row.names(records) <- NULL

# Sum of mutations per patient
records$sum_mutations <- rowSums(records[,c("filtered_SNVs_common", "filtered_indels_common")])

# Removing carcinoid histology patients
records <- records[!records$histology == 4,]
records <- records[!records$Patient_ID == "NSLC-0124",]

# Dataset for the histology plot
histology_mutations <-aggregate(records$sum_mutations, by=list(records$histology), FUN=sum)
colnames(histology_mutations)<- c("histology_type", "total_mutations")
histology_mutations$histology_type[histology_mutations$histology_type == 2] <- "squamous cell carcinoma (2)"
histology_mutations$histology_type[histology_mutations$histology_type == 3] <- "adenocarcinoma (3)"
#histology_mutations$histology_type[histology_mutations$histology_type == 4] <- "carcinoid tumor (4)"
histology_mutations$histology_type[histology_mutations$histology_type == 6] <- "adenosquamous carcinoma (6)"
histology_mutations$histology_type[histology_mutations$histology_type == 7] <- "sarcomatoid carcinoma (7)"

# Making of plots
Palette <- c("#a63a9b","#3D87B8","#2a9d8f","#9c818e","#e9c46a","#dc8f41","#Fe5bac","#F50025")

# Dataset for the patient plot
patient_df <- reshape2::melt(records[,c("Patient_ID", "coding_mutations", "noncoding_mutations")], id.vars="Patient_ID", variable.name = "Nature", value.name = "Count")

# Patient plot divided by coding and non-coding mutations
ggplot(patient_df, aes(x=reorder(substring(patient_df[,"Patient_ID"], 6, 10), -Count), y=Count, fill=Nature)) +
  #ggtitle("Patients' filtered mutation count") +
  geom_bar(stat = "identity", width=0.9, position = "stack") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5, colour = "#264653"), axis.text.y = element_text(colour = "#264653"), axis.ticks.x=element_blank(),
        panel.grid.major.y = element_line(colour="black", size=0.5), legend.position = c(0.7,0.65), legend.text=element_text(size=11), plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  scale_fill_manual(values = Palette[c(3,5)], labels = c("Coding", "Non-coding")) +
  scale_y_continuous(name="Number of variants", expand = c(0, 0)) +
  scale_x_discrete(name="Patients")

ggsave("nb_variants_patients.png", width = 12, height = 6, dpi = 1200)

# Dataset for the histology plot with Number of variants divided in coding and non-coding mutations
histology_df <- reshape2::melt(records[,c("Patient_ID", "coding_mutations", "noncoding_mutations", "histology")], id.vars=c("Patient_ID", "histology"), variable.name = "Nature", value.name = "Count")
histology_df$histology = factor(histology_df$histology, levels=c('3','2','4','6','7'))
# Corresponding levels : 3 = adenocarcinoma, 2 = squamous cell carcinoma, 4 = carcinoid tumor, 6 = adenosquamous carcinoma, 7 = sarcomatoid carcinoma
levels(histology_df$histology) <- c('adenocarcinoma', 'squamous cell\ncarcinoma', 'carcinoid\ntumor', 'adenosquamous\ncarcinoma', 'sarcomatoid\ncarcinoma')

# Histology plot with Number of variants divided in coding and non-coding mutations
reorder(substring(histology_df[which(histology_df$histology=='squamous cell\ncarcinoma'),"Patient_ID"], 6, 10),
        -histology_df$Count[which(histology_df$histology=='squamous cell\ncarcinoma')]) # Patients 20 and 139 were not shown in the right order.

png("histology_v2.png", width = 18, height = 9, units = "in", res = 1200)
# There is a bad interaction between geom_col and log scales (wanted for y-axis). So instead of doing a simple transformation of the y-axis inside "scale_y_continuous()",
# a y-axis break is done to show coding mutations (from 0 to 500 mutations) and non coding mutations (from 0 to 40 000 mutations).
{
  top <- ggplot(histology_df, aes(x=reorder(substring(histology_df[,"Patient_ID"], 6, 10), -Count), y=Count, fill=Nature)) +
    #ggtitle("Distribution of histological types' filtered mutation count") +
    geom_col(width=0.9, position = position_stack(reverse = TRUE)) + 
    theme_classic() +
    facet_grid(~histology, scales = "free_x", space = "free_x", switch = "x") +
    theme(axis.text.y = element_text(colour = "#264653", size = 12),
          panel.grid.major.y = element_line(colour="black", size=0.5), legend.position = c(0.8,0.47), legend.text=element_text(size=16), plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(), panel.spacing = unit(0.5, "lines"), plot.margin=unit(c(5.5,5.5,0,5.5),"pt"), axis.title.y = element_text(hjust = 0, size = 18),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          strip.text.x=element_blank()) +
    scale_fill_manual(values = Palette[c(3,5)], labels = c("Coding", "Non-coding"), guide = guide_legend(reverse = TRUE)) +
    scale_color_manual(values = Palette) +
    scale_y_continuous(name="Number of variants", expand = c(0, 0)) +
    coord_cartesian(ylim = c(500,NA), clip = "on")
  
  top.gp <- ggplotGrob(top)
  
  # Bottom zooming on coding mutations and showing the x-axis with histology names.
  bot <- ggplot(histology_df, aes(x=reorder(substring(histology_df[,"Patient_ID"], 6, 10), -Count), y=Count, fill=Nature)) +
    #ggtitle("Distribution of histological types' filtered mutation count") +
    geom_col(width=0.9, position = position_stack(reverse = TRUE)) + 
    theme_classic() +
    facet_grid(~histology, scales = "free_x", space = "free_x", switch = "x") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, colour = "#264653"), axis.text.y = element_text(colour = "#264653", size = 12), axis.ticks.x=element_blank(),
          panel.grid.major.y = element_line(colour="black", size=0.5), legend.position = "none", legend.text=element_blank(), legend.background = element_blank(), 
          plot.title = element_text(hjust = 0.5), legend.title=element_blank(), panel.spacing = unit(0.5, "lines"), strip.placement = "outside", 
          strip.text.x = element_text(angle = 90, colour = "#264653", size = 14), axis.title.y = element_blank(), axis.title.x = element_text(size = 18),
          strip.switch.pad.grid = unit(-1, "cm"), plot.margin=unit(c(5,5.5,45,5.5),"pt")) +
    scale_fill_manual(values = Palette[c(3,5)]) +
    scale_color_manual(values = Palette) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(name="Patients") + 
    coord_cartesian(ylim = c(0,500), clip = "on")
  
  bot.gp <- ggplotGrob(bot)
  
  t_title <- bot.gp$layout[bot.gp$layout[['name']] == 'xlab-b' ,][['b']]
  t_strip <- bot.gp$layout[grepl('strip', bot.gp$layout[['name']]),][['b']]
  
  bot.gp$layout[bot.gp$layout[['name']] == 'xlab-b' ,][['t']] <- unique(t_strip)
  bot.gp$layout[bot.gp$layout[['name']] == 'xlab-b' ,][['b']] <- unique(t_strip)
  bot.gp$layout[grepl('strip', bot.gp$layout[['name']]),][['t']] <- t_title
  bot.gp$layout[grepl('strip', bot.gp$layout[['name']]),][['b']] <- t_title
  
  for(i in which(grepl("strip-b", bot.gp$layout$name))){
    bot.gp$grobs[[i]]$layout$clip <- "off"
  }
  
  for (k in grep("strip-b",bot.gp$layout$name)) {
    bot.gp$grobs[[k]]$grobs[[1]]$children[[1]] <- grid.lines(x=unit(c(0,1),"npc"), y=unit(c(1,1),"npc"), 
                                                             gp=gpar(col="#2a9d8f", lwd=2))
  }
  
  # Assembly of the bottom and top plots.
  plot_grid(top.gp, bot.gp, rel_heights=c(0.5, 0.5), ncol = 1, align = "v", vjust = 0)
}
dev.off()

# Previous version of the histology plot (without zooming on coding mutations)
png("histology.png", width = 12, height = 6, units = "in", res = 1200)
{ histology_plot <- ggplot(histology_df, aes(x=reorder(substring(histology_df[,"Patient_ID"], 6, 10), -Count), y=Count, fill=Nature)) +
    #ggtitle("Distribution of histological types' filtered mutation count") +
    geom_col(width=0.9, position = position_stack(reverse = FALSE)) + # reverse = TRUE to have coding mutations bars under non-coding mutations bars.
    theme_classic() +
    facet_grid(~histology, scales = "free_x", space = "free_x", switch = "x") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6, colour = "#264653"), axis.text.y = element_text(colour = "#264653"), axis.ticks.x=element_blank(),
          panel.grid.major.y = element_line(colour="black", size=0.5), legend.position = c(0.8,0.65), legend.text=element_text(size=11), plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(), panel.spacing = unit(0.5, "lines"), strip.placement = "outside", strip.text.x = element_text(angle = 90, colour = "#264653", size = 10), 
          strip.switch.pad.grid = unit(-1, "cm"), plot.margin=unit(c(5.5,5.5,35,5.5),"pt")) +
    scale_fill_manual(values = Palette[c(3,5)], labels = c("Coding", "Non-coding"), guide = guide_legend(reverse = TRUE)) +
    scale_color_manual(values = Palette) +
    scale_y_continuous(name="Number of variants", expand = c(0, 0)) +
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


# Dataset for the histology plot with Number of variants divided in SNPs and indels
histology_df.2 <- reshape2::melt(records[,c("Patient_ID", "filtered_indels_common", "filtered_SNVs_common", "histology")], id.vars=c("Patient_ID", "histology"), variable.name = "Nature", value.name = "Count")
histology_df.2$histology = factor(histology_df.2$histology, levels=c('3','2','4','6','7'))
# Corresponding levels : 3 = adenocarcinoma, 2 = squamous cell carcinoma, 4 = carcinoid tumor, 6 = adenosquamous carcinoma, 7 = sarcomatoid carcinoma
levels(histology_df.2$histology) <- c('adenocarcinoma', 'squamous cell\ncarcinoma', 'carcinoid\ntumor', 'adenosquamous\ncarcinoma', 'sarcomatoid\ncarcinoma')

# Histology plot with number of mutation divided in SNPs and indels
reorder(substring(histology_df.2[which(histology_df.2$histology=='squamous cell\ncarcinoma'),"Patient_ID"], 6, 10),
        -histology_df.2$Count[which(histology_df.2$histology=='squamous cell\ncarcinoma')]) # Patients 20 and 139 were not shown in the right order.
png("histology_SNPs-indels.png", width = 12, height = 6, units = "in", res = 1200)
{histology_plot <- ggplot(histology_df.2, aes(x=reorder(substring(histology_df.2[,"Patient_ID"], 6, 10), -Count), y=Count, fill=Nature)) +
    #ggtitle("Distribution of histological types' filtered mutation count") +
    geom_col(width=0.9, position = "stack") +
    theme_classic() +
    facet_grid(~histology, scales = "free_x", space = "free_x", switch = "x") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6, colour = "#264653"), axis.text.y = element_text(colour = "#264653"), axis.ticks.x=element_blank(),
          panel.grid.major.y = element_line(colour="black", size=0.5), legend.position = c(0.8,0.65), legend.text=element_text(size=11), plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(), panel.spacing = unit(0.5, "lines"), strip.placement = "outside", strip.text.x = element_text(angle = 90, colour = "#264653", size = 10), 
          strip.switch.pad.grid = unit(-1, "cm"), plot.margin=unit(c(5.5,5.5,35,5.5),"pt")) +
    scale_fill_manual(values = Palette[c(3,5)], labels = c("Indels", "SNPs")) +
    scale_color_manual(values = Palette) +
    scale_y_continuous(name="Number of variants", expand = c(0, 0)) +
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

