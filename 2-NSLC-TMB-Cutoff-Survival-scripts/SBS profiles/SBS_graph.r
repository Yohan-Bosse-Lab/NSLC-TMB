library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(viridis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

SBS <- read.table("NSLC_SBS/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt", sep = "\t", header = TRUE)
SBS.df <- reshape2::melt(SBS, id.vars="Samples", variable.name = "SBS_type", value.name = "Count")
SBS.df$Samples<-gsub("NSLC_","",as.character(SBS.df$Samples))

ggplot(SBS.df, aes(x=reorder(Samples, -Count), y=Count, fill=SBS_type)) +
  geom_bar(position="stack", stat="identity", width = 0.9) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6), axis.ticks.x=element_blank(),
        panel.grid.major.y = element_line(colour="black", linewidth=0.5), legend.text=element_text(size=11), axis.title = element_text(size=12),
        legend.title=element_blank(), panel.spacing = unit(0.5, "lines"), strip.placement = "outside", strip.text.x = element_text(angle = 90, size = 10), 
        plot.margin=unit(c(5.5,5.5,35,5.5),"pt"), panel.grid.minor = element_blank()) +
  scale_y_continuous(name="Number of mutations in each signature",  limits = c(0, NA), expand = c(0, 0)) +
  scale_x_discrete(name="Patients") +
  #scale_fill_viridis(discrete = TRUE, option = "D")
  scale_fill_ucscgb()
  #scale_fill_brewer(palette = "RdYlBu")
  
ggsave(file="COSMIC_SBS96.png", width=16, height=8, dpi=300)
  