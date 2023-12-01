
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(survminer)
library(DescTools)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################## Data import and formatting ############################
################################################################################

## Difference in TMB between patients with and without event within 2 years
# Adjusting PFS event : death or cancer progression within 2 years of surgery
studied_patients_data <- fread("Results/Feature Selection/studied_patients_data_5k_50k_500k_TMB.csv")
nrow(studied_patients_data[PFS_indicator_adjusted == 1])
nrow(studied_patients_data[PFS_indicator_adjusted != 1])

# Fetching and merging all patients variants
VEP_data_all_patients <- fread("Data/n82_NSLC_variants_data.csv", keepLeadingZeros = TRUE)

## Kaplan-Meier plt of OS and PFS with dicotomous TMB
OS.dicho <- survfit(coxph(Surv(time=OS_time, VitalStatus) ~ strata(TMB_high_low) + pathological_stage_refactor, data = studied_patients_data))
PFS.dicho <- survfit(coxph(Surv(time=PFS_time, PFS_indicator) ~ strata(TMB_high_low) + pathological_stage_refactor, data = studied_patients_data))

OS_PFS_surv <- ggsurvplot_combine(list(OS.dicho, PFS.dicho), 
                                  ggtheme=theme_minimal(), 
                                  data = clin_data, 
                                  palette = c("#457B9D","#1D3557", "#EC9A9A", "#E63946"), 
                                  risk.table = T,
                                  xlab = "Time after surgery (days)",
                                  ylab = "Proportion surviving",
                                  break.x.by = 500,
                                  surv.median.line = "hv",
                                  legend.title = "",
                                  legend.labs = c("OS High TMB", "OS Low TMB", "PFS High TMB", "PFS Low TMB"))

OS_PFS_surv$plot <- OS_PFS_surv$plot + geom_vline(xintercept=1826, linetype = "solid") + # 5 years
  geom_vline(xintercept=730, linetype = "solid") + # 2 years
  geom_segment(x=5466, y = 0.5, yend = 0, xend = 5466, linetype = "dashed") + # OS High TMB median
  annotate("text", x = 5466, y=0.005, label = "5466", cex=3.3, vjust=0, hjust = 1.1) + # OS High TMB median label
  geom_segment(x=2829, y = 0.5, yend = 0, xend = 2829, linetype = "dashed") + # PFS High TMB median
  annotate("text", x = 2829, y=0.005, label = "2829", cex=3.3, vjust=0, hjust = 1.1) + # PSS High TMB median label
  geom_segment(x=4147, y = 0.5, yend = 0, xend = 4147, linetype = "dashed") + # PFS Low TMB median
  annotate("text", x = 4147, y=0.005, label = "4147", cex=3.3, vjust=0, hjust = 1.1) # PSS Low TMB median label

png("Poster FMED/OS_PFS.png", width = 8, height = 8, units = "in", res = 1200)
print(OS_PFS_surv)
dev.off()

# Median of events for different genome regions
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, complete_WGS_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, complete_WES_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, TMB_per_region_exons_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, TMB_per_region_introns_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, TMB_per_region_intergenic_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, TMB_per_region_regulatory_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, sw_TMB_cv_50k])
MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, sw_TMB_cv_5k])

# Median of non-events for different genome regions
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, complete_WGS_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, complete_WES_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, TMB_per_region_exons_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, TMB_per_region_introns_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, TMB_per_region_intergenic_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, TMB_per_region_regulatory_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, sw_TMB_cv_50k])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, sw_TMB_cv_5k])

# Mann-Whitney U Test for the different genome regions
wilcox.test(complete_WGS_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(complete_WES_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(TMB_per_region_exons_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(TMB_per_region_introns_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(TMB_per_region_intergenic_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(TMB_per_region_regulatory_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(sw_TMB_cv_50k ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value
wilcox.test(sw_TMB_cv_5k ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value

################################################################################
########################## COSMIC LUAD genes TMB ###############################
################################################################################
library(mygene)
library(AnnotationHub)

# Genes acquired from https://cancer.sanger.ac.uk/cosmic/browse/tissue?wgs=off&sn=lung&ss=all&hn=carcinoma&sh=non_small_cell_carcinoma&in=t&src=tissue&all_data=n
COSMIC_genes.names <- c("EGFR", "KRAS", "TP53", "ARID1A","PIK3CA","ALK","BRAF",
                        "MET","STK11","KEAP1","CDKN2A","ATM","KMT2D","SMARCA4",
                        "NF1","SETBP1","FAT1","ERBB2","TET2","NOTCH1","NTRK1", 
                        "NTRK2", "NTRK3","RET","ROS1")

# Acquire Ensembl IDs
COSMIC_genes.ensembl_id <- unique(na.omit(queryMany(COSMIC_genes.names, scopes = "symbol", species = "human",  fields="ensembl.gene"))$ensembl.gene)

# Build GRCh37 database
ah = AnnotationHub()
#AnnotationHub::query(ah, c('gtf', 'Homo_sapiens', 'GRCh37'))
GRCh37.gtf <- ah[['AH10684']]

## subset the gtf files for only protein_coding genes and lincRNAs
GRCh37.gtf<- GRCh37.gtf[GRCh37.gtf$gene_biotype %in% c('protein_coding', 'lincRNA')]
table(GRCh37.gtf$gene_biotype)

## make a txdb and keep conventional chromosomes
GRCh37.txdb <- makeTxDbFromGRanges(GRCh37.gtf)
GRCh37.txdb <- keepSeqlevels(GRCh37.txdb, c(1:22, 'X', 'Y'), pruning.mode = 'coarse')

## cds (not compiled in the total genome size because included in exons)
cds <- cdsBy(GRCh37.txdb, 'gene') %>% unstrand() %>% GenomicRanges::reduce()

# Get genes CDS size with GRCh37
COSMIC_genes.size <- setNames(sapply(COSMIC_genes.ensembl_id, function(x) sum(width(cds[[x]]))), nm=COSMIC_genes.names)

# Table summarizing the ensembl id, the gene name and its size
COSMIC_genes <- data.table(ensembl_id = COSMIC_genes.ensembl_id,
                           names = COSMIC_genes.names,
                           size = COSMIC_genes.size)



# WGS data with annotation for the genes of interest
COSMIC_genes.data <- setNames(lapply(COSMIC_genes.names, function(gene) VEP_data_all_patients[SYMBOL %like% gene]), nm = COSMIC_genes.names) %>% rbindlist()

# Patients list
Patients_list <- gsub("NSLC-", "",sort(studied_patients_data[, Patient_ID]))

# Mutation count for each patient of the genes of interest
studied_patients_data[, COSMIC_mutation_count := sapply(Patients_list, function(patient) nrow(COSMIC_genes.data[as.numeric(Patient_ID) == patient]))]
studied_patients_data[, COSMIC_TMB := COSMIC_mutation_count / (sum(COSMIC_genes$size)/1e6)]

MedianCI(studied_patients_data[PFS_indicator_adjusted == 1, COSMIC_TMB])
MedianCI(studied_patients_data[PFS_indicator_adjusted != 1, COSMIC_TMB])
wilcox.test(COSMIC_TMB ~ PFS_indicator_adjusted, data=studied_patients_data, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)$p.value

fwrite(studied_patients_data, file = "Results/Feature Selection/studied_patients_data_5k_50k_500k_COSMIC_TMB.csv")

################################################################################
################################ PFS graphs ####################################
################################################################################

ggplot(clin_data[!(PFS_indicator == 0 & pfs_two == "<2")], aes(x=pfs_two, y=complete_WGS_TMB)) +
  geom_boxplot(data = clin_data[PFS_indicator == 1 & pfs_two == "<2"], width=0.5, outlier.colour = "white") +
  geom_boxplot(data = clin_data[pfs_two == ">=2"], width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor),
              size = 2, width = 0.15) +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  ylab("TMB g?nome entier (mut/Mb)") + 
  scale_x_discrete(labels = c("<2" = "<2 ans\nn=10",
                              ">=2" = "\u22652 ans\nn=63")) +
  scale_y_continuous(breaks = seq(0,10,1),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_manual(name = "Pathological stage", values = c("#1D3557", "#EC9A9A", "#E63946")) +
  # annotate("text", x=0.5, y = 1.55, label = "n=8", hjust = "left") +
  # annotate("text", x=0.5, y = 1.90, label = "n=12", hjust = "left") +
  # annotate("text", x=2.4, y = 1.55, label = "n=42", hjust = "left") +
  # annotate("text", x=2.4, y = 1.90, label = "n=18", hjust = "left") +
  # annotate("text", x=0, y = 10, label = glue("Mann-Whitney p={WGS.two.WRS.p}"), hjust = -0.05) +
  # annotate("text", x=0, y = 9.5, label = glue("Event Mann-Whitney p={WGS.two.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off") +
  theme(legend.position=c(.2,.85),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=11),
        axis.title.y = element_text(size=13),
        axis.text = element_text(size=13, color = "black"),
        axis.title.x = element_blank())

ggsave(file=glue("Poster FMED/PFS_WGS_2years_event.png"), width=5, height=5, dpi=300)


################################################################################
############################ Median TMB Violon plot ############################
################################################################################

# Importing the dataset containing all TMB calculations.
studied_patients_data <- fread("Results/Feature Selection/studied_patients_data_5k_50k_500k_COSMIC_TMB.csv")

# Transforming the data for violon plot purposes
studied_patients_data.melt <- melt(studied_patients_data[, .(PFS_indicator_adjusted, 
                                                             complete_WGS_TMB, 
                                                             TMB_per_region_intergenic_TMB, 
                                                             TMB_per_region_introns_TMB,
                                                             TMB_per_region_regulatory_TMB,
                                                             TMB_per_region_exons_TMB,
                                                             sw_TMB_cv_50k,
                                                             sw_TMB_cv_5k,
                                                             COSMIC_TMB)], 
                                   id.vars = "PFS_indicator_adjusted", variable.name = "Genome region", value.name = "TMB")

# Changing PFS_indicator_adjusted to be a factor instead of a numerical binary variable
studied_patients_data.melt[, PFS_indicator_adjusted := as.factor(PFS_indicator_adjusted)]

# Special violon plot where the groups' plot are splited and put next to each other.
# Source: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(self,
                        data,
                        ...,
                        # add the nudge here
                        nudge = 0,
                        draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- plyr::arrange(transform(data,
                                       x = if (grp %% 2 == 1) xminv else xmaxv),
                             if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ],
                     newdata,
                     newdata[nrow(newdata), ],
                     newdata[1, ])
    newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    
    # now nudge them apart
    newdata$x <- ifelse(newdata$group %% 2 == 1,
                        newdata$x - nudge,
                        newdata$x + nudge)
    
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      
      quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                           draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)),
                         setdiff(names(data), c("x", "y")),
                         drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                      quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       ggplot2::GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              # nudge param here
                              nudge = 0,
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  
  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = stat,
                 geom = GeomSplitViolin,
                 position = position,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(trim = trim,
                               scale = scale,
                               # don't forget the nudge
                               nudge = nudge,
                               draw_quantiles = draw_quantiles,
                               na.rm = na.rm,
                               ...))
}

# Region names for the violon plot
genome_region_names <- c(
  "complete_WGS_TMB" = "Whole genome\n(~3100 Mb)",
  "TMB_per_region_intergenic_TMB" = "Intergenic\n(~1500 Mb)",
  "TMB_per_region_introns_TMB" = "Introns\n(~1300 Mb)",
  "TMB_per_region_regulatory_TMB" = "Regulatory\n(~230 Mb)",
  "TMB_per_region_exons_TMB"  = "Whole exome\n(~35 Mb)",
  "sw_TMB_cv_50k" = "Optimized 0.05 Mb\n(~20.73 Mb)",
  "sw_TMB_cv_5k" = "Optimized 0.005 Mb\n(~3.94 Mb)",
  "COSMIC_TMB" = "NSCLC genes\n(~0.13 Mb)")

## P-value annotation positions
## Source: https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
# Getting the y position where to put the annotation for each region. This value is represented by the max TMB value for each region.
max_regions_value <- sapply(unique(studied_patients_data.melt$`Genome region`), 
                            function(region) max(studied_patients_data.melt[`Genome region` == region, TMB]), 
                            USE.NAMES = TRUE)

dat_text <- data.table(label = c("p-value = 0.015*", "p-value = 0.013*", "p-value = 0.016*", "p-value = 0.020*", "p-value = 0.107", "p-value < 0.001***", "p-value < 0.001***", "p-value = 0.070"),
                       "Genome region" = unique(studied_patients_data.melt$`Genome region`),
                       y = max_regions_value,
                       PFS_indicator_adjusted = factor(rep(0,8)))

# Making the violon plot, fill = PFS_indicator_adjuste
ggplot(data=studied_patients_data.melt, aes(x=`Genome region`, y=TMB, fill = PFS_indicator_adjusted)) +
  theme_classic() +
  geom_split_violin(nudge =.035, width = 0.9, linewidth = 0) +
  geom_half_boxplot(data=studied_patients_data.melt[PFS_indicator_adjusted==1],
                    fill = "white", width = 0.15, color = "#103F59", 
                    outlier.colour = "#103F59", outlier.size = 1, size = 0.4,
                    position = position_nudge(x=.035), side = "r", errorbar.draw = TRUE, errorbar.length = 1) +
  geom_half_boxplot(data=studied_patients_data.melt[PFS_indicator_adjusted==0], 
                    fill = "white", width = 0.15, color = "#E9C46A", 
                    outlier.colour = "#E9C46A", outlier.size = 1, size = 0.4,
                    position = position_nudge(x=-.035), errorbar.draw = TRUE, errorbar.length = 1) +
  scale_y_continuous("TMB, log10-scaling", 
                     breaks=c(0, .25, .5, 1, 2.5, 5, 10, 25, 50, 100, 250),
                     limits = c(0,430),
                     trans=scales::pseudo_log_trans(base = 10),
                     expand = c(0,0)) +
  scale_fill_manual(name = "Patients", 
                    labels = c("Event, n=10",
                               "No event, n=63"), 
                    values = c("#103F59", "#E9C46A")) +
  facet_grid(~`Genome region`, 
             scales = "free_x", 
             space = "free_x", 
             switch = "x", 
             labeller = labeller(.default = capitalize,`Genome region` = genome_region_names)) +
  geom_text(size = 3.5,
            data = dat_text,
            mapping = aes(x = 1, y = y, label = label),
            hjust = 0.5,
            vjust = -1.05) +
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),  # Make facet label background white.
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/Model comparison/genome_regions.png"), width=12, height=6, dpi=500)
