
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(karyoploteR)
library(regioneR)
library(ggbio)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/n82_NSLC_surv_data.xlsx", sheetIndex = 1)) # See script Data_harmonization.r for details

# Fetching all patients variants
VEP_data_all_patients <- fread("Data/n82_NSLC_variants_data.csv") # See script Data_harmonization.r for details

################################################################################
############  CLASSIC TMB METHODS TO BE OPTIMIZED: WGS AND WES #################
################################################################################

# WGS TMB mean and median for patients with event (progression or death) before 2 years following surgery
mean(clin_data[PFS_indicator == 1 & pfs_two == "<2", complete_WGS_TMB])
median(clin_data[PFS_indicator == 1 & pfs_two == "<2", complete_WGS_TMB])

# WGS TMB mean and median for the other patients (event or not after two years)
mean(clin_data[pfs_two != "<2", complete_WGS_TMB])
median(clin_data[pfs_two != "<2", complete_WGS_TMB])

# WES TMB mean and median for patients with event (progression or death) before 2 years following surgery
mean(clin_data[PFS_indicator == 1 & pfs_two == "<2", complete_WES_TMB])
median(clin_data[PFS_indicator == 1 & pfs_two == "<2", complete_WES_TMB])

# WES TMB mean and median for the other patients (event or not after two years)
mean(clin_data[pfs_two != "<2", complete_WES_TMB])
median(clin_data[pfs_two != "<2", complete_WES_TMB])

################################################################################
################ USING CADD-PHRED SCORE TO CHOOSE HOTSPOTS #####################
################################################################################

## Datasets for the analyses
# Patients with event: progression or death before 2 years following surgery
event_patients.list <- gsub("NSLC-", "", clin_data[PFS_indicator == 1 & pfs_two == "<2", Patient_ID])
event_patients.variants <- rbindlist(lapply(event_patients.list, function(x) VEP_data_all_patients[str_detect(VEP_data_all_patients[,ID], glue("^{x}:")),]), use.names = TRUE)
event_patients.clin <- rbindlist(lapply(event_patients.list, function(x) clin_data[str_detect(clin_data[,Patient_ID], glue("NSLC-{x}")),]), use.names = TRUE)

# Patients without event: no progression or no death before 2 years following surgery
no_event_patients.list <- gsub("NSLC-", "", clin_data[pfs_two != "<2", Patient_ID])
no_event_patients.variants <- rbindlist(lapply(no_event_patients.list, function(x) VEP_data_all_patients[str_detect(VEP_data_all_patients[,ID], glue("^{x}:")),]), use.names = TRUE)
no_event_patients.clin <- rbindlist(lapply(no_event_patients.list, function(x) clin_data[str_detect(clin_data[,Patient_ID], glue("NSLC-{x}")),]), use.names = TRUE)

## Vizualisation
# Top 10% PHRED scaled CADD score cutoff for hotspot regions
hotspot_cutoff.CADD <- 20

# Arranging chromosome names for ggplot
chrom_name <- c(1:22, 'X', 'Y')
event_patients.variants$CHROM <- factor(event_patients.variants$CHROM, levels=chrom_name)
no_event_patients.variants$CHROM <- factor(no_event_patients.variants$CHROM, levels=chrom_name)

# Plot for positive event
ggplot(event_patients.variants,
       aes(x=POS, y=CADD_PHRED, color = CHROM)) +
  geom_line(linewidth = 0.2) +
  geom_hline(yintercept = hotspot_cutoff.CADD, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside', # for facet_grid label position
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none") + # for facet_grid space between panels
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,5)) +
  xlab("Chromosomes") + 
  ylab("PHRED CADD score") +
  facet_grid(~CHROM,
             space = "free_x",
             scales = "free_x",
             switch = "x")

ggsave(file=glue("Results/CADD hotspots/CADD_score_event.png"), width=15, height=7, dpi=300)

# Plot for negative event
ggplot(no_event_patients.variants, 
       aes(x=POS, y=CADD_PHRED, color = CHROM)) +
  geom_line(linewidth = 0.2) +
  geom_hline(yintercept = hotspot_cutoff.CADD, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside', # for facet_grid label position
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none") + # for facet_grid space between panels
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,5)) +
  xlab("Chromosomes") + 
  ylab("PHRED CADD score") +
  facet_grid(~factor(CHROM, levels=chrom_name),
             space = "free_x",
             scales = "free_x",
             switch = "x")

ggsave(file=glue("Results/CADD hotspots/CADD_score_no_event.png"), width=15, height=7, dpi=300)


## CHOICE OF HOTSPOTS
high_CADD_mut  <- event_patients.variants[CADD_PHRED >= hotspot_cutoff.CADD,]
region_bounds <- 100 # amount of bp to look around the high CADD scoring mutation

high_CADD_mut.regions <- rbindlist(lapply(1:nrow(high_CADD_mut), function(x) event_patients.variants[POS >= high_CADD_mut[x, POS - region_bounds] & END_POS <= high_CADD_mut[x, END_POS + region_bounds], ]))
high_CADD_mut.regions.size  <- sapply(1:nrow(high_CADD_mut), function(x) high_CADD_mut[x, END_POS + region_bounds] - high_CADD_mut[x, POS - region_bounds])

total_size <- sum(high_CADD_mut.regions.size)

nrow(high_CADD_mut.regions) / total_size # The final TMB for the regions around the high scoring CADD mutations

summary(coxph(data = clin_data, Surv(time=time_RpFS, event_relapse_two_years_indicator) ~ TMB_high_low))$waldtest[3]
colnames(clin_data)

################################################################################
################ USING VARIANT FREQUENCY TO CHOOSE HOTSPOTS ####################
################################################################################

## TODO
# Use specific highly mutated regions of positive patients (event_relapse_two_years_indicator == 1) to optimize rather than overall highly mutated regions

# Window size used to compute sliding windows
window_size <- 5000
Mb_length <-  window_size/1e6

## Merging all chromosome sliding windows densities together
chrom_name <- c(1:22, 'X', 'Y')
# Reading density data of chromosomes with a specific window size
chrom_density <- rbindlist(lapply(chrom_name, 
                                  function(x) fread(glue("Data/Sliding windows and frequency/chromosome_{x}_sw_density_{Mb_length}Mb.csv"))[, chrom:=glue("{x}")]),
                           use.names = TRUE)

# Adjusting the class required for manipulation
attr(chrom_density$chrom, 'glue') <- NULL

# *** Percentage of highest variant frequency windows to keep
top_windows_percentage <- 0.02 # Keep top 0.1% of windows with the most variants

# Get the windows with the most variants
setorder(chrom_density, -variant_density)
hotspot.cutoff <- ceiling(nrow(chrom_density) * top_windows_percentage)
hotspot.windows <- chrom_density[1:hotspot.cutoff]

# Merge the neighboring windows (combined neighboring windows)
hotspot.windows.reduced <- GRanges(hotspot.windows) %>% GenomicRanges::reduce() %>% as.data.table() %>% .[, strand:=NULL]
names(hotspot.windows.reduced) <- c("CHROM", "POS", "END_POS", "WIDTH")


### Merged neighboring windows analysis
# Extract the variants inside each merged window
hotspots.list.merge <- lapply(1:nrow(hotspot.windows.reduced), 
                              function(x){rbindlist(list(hotspot.windows.reduced[x], 
                                                         setorder(VEP_data_all_patients[POS >= hotspot.windows.reduced[x, POS] & 
                                                                                          END_POS <= hotspot.windows.reduced[x, END_POS] & 
                                                                                          CHROM == hotspot.windows.reduced[x, CHROM]], POS)),
                                                    fill=TRUE)
                              }
)

# Total merged windows size, needed for TMB calculation
hotspots.list.merged.total_size <- sum(hotspot.windows.reduced[, WIDTH])

# Order the merged windows based on the number of variants they include
hotspots.list.merged.nrow <- sapply(1:length(hotspots.list.merge), function(x) nrow(hotspots.list.merge[[x]]))
hotspots.list.merged.ordered <- hotspots.list.merge[order(hotspots.list.merged.nrow, decreasing = TRUE)]

# Fetching the number of variants included in the merged windows hotspots for each patient
hotspots.list.merged.variant_count.list <- setNames(lapply(gsub("NSLC-", "", clin_data[, Patient_ID]), 
                                                           function(patient) sum(unlist(lapply(1:length(hotspots.list.merged.ordered),
                                                                                               function(hotspot) sum(na.omit(str_detect(hotspots.list.merged.ordered[[hotspot]][, ID], glue("^{patient}:")))))))),
                                                    nm = clin_data[, Patient_ID])

# Final variant count and TMB for each patients based on the mutation frequency hotspots
hotspots.list.merged.variant_count <- data.table(Patient_ID = names(hotspots.list.merged.variant_count.list),
                                                 Merged_window.Variant_count = unname(unlist(hotspots.list.merged.variant_count.list)),
                                                 Merged_window.TMB = unname(unlist(hotspots.list.merged.variant_count.list)) / hotspots.list.merged.total_size)

# Integrating hotspot variants count and TMB with the rest of the data
clin_data.variant_count <- merge(clin_data, hotspots.list.merged.variant_count, by.x = "Patient_ID", by.y = "Patient_ID")

# Temporary. Trying to maximize the LRT p-value on the chosen windows
library(survival)
colnames(clin_data.variant_count)
summary(coxph(data = clin_data.variant_count, Surv(time=time_RpFS, event_relapse_two_years_indicator) ~ TMB_high_low))$waldtest[3]
coxph(data = clin_data.variant_count, Surv(time=time_RpFS, event_relapse_two_years_indicator) ~ pathological_stage_refactor)
summary(coxph(data = clin_data.variant_count, Surv(time=time_RpFS, event_relapse_two_years_indicator) ~ complete_WGS_TMB))$waldtest[3]
summary(coxph(data = clin_data.variant_count, Surv(time=time_RpFS, event_relapse_two_years_indicator) ~ Merged_window.Variant_count))$waldtest[3]

### Unmerged neighboring windows analysis
# If windows are not be merged, change the names of the columns to avoid names matching problems later on
names(hotspot.windows) <- c("POS", "END_POS", "FREQUENCY", "CHROM")
hotspot.windows[, CHROM:=as.character(CHROM)] # Editing 

# Extract the variants inside each window
hotspots.list <- lapply(1:nrow(hotspot.windows), 
                        function(x){rbindlist(list(hotspot.windows[x], 
                                                   setorder(VEP_data_all_patients[POS >= hotspot.windows[x, POS] & 
                                                                                    END_POS <= hotspot.windows[x, END_POS] & 
                                                                                    CHROM == hotspot.windows[x, CHROM]], POS)),
                                              fill=TRUE)
                        }
)

# Order the windows based on the number of variants they include
hotspots.list.nrow <- sapply(1:length(hotspots.list), function(x) nrow(hotspots.list[[x]]))
hotspots.list.ordered <- hotspots.list[order(hotspots.list.nrow, decreasing = TRUE)]



# Looking for specific genes
min(VEP_data_all_patients[str_detect(VEP_data_all_patients[, SYMBOL], "EGFR"), POS])
max(VEP_data_all_patients[str_detect(VEP_data_all_patients[, SYMBOL], "EGFR"), END_POS])
unique(str_extract(VEP_data_all_patients[str_detect(VEP_data_all_patients[, SYMBOL], "EGFR"), ID],"[^:]+"))



################################################################################
################################ HOTSPOT PLOTS #################################
################################################################################

# Reference: https://stackoverflow.com/questions/58696329/is-it-possible-to-over-ride-the-x-axis-range-in-r-package-ggbio-when-using-autop
library(ggbio)
library(ggplotify)
library(gridExtra)
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75

# Some functions for plot x axis ticks
nearest_lower_limit <- function(value, divider) {
  value - (value %% divider)
}
nearest_upper_limit <- function(value, divider) {
  value + (value %% divider)
}


### Unmerged neighboring windows
# Specific variant chosen for first test. 
mut.7.142465001 <- hotspots.list.ordered[1] %>% makeGRangesFromDataFrame(., ignore.strand = TRUE,
                                                                         start.field="POS",
                                                                         end.field="END_POS",
                                                                         strand.field="CHROM")
mut.7.142465001.dt <- as.data.table(mut.7.142465001)

freq <- ggplot(data = mut.7.142465001.dt[c(2:.N)], aes(x=start)) +
  #geom_point(data = as.data.table(mut.7.142465001)[2:length(mut.7.142465001)], aes(x=start, y=1)) +
  geom_histogram(binwidth = 1) +
  theme_classic()

endb.anno <- autoplot(ensdb, which = mut.7.142465001) + theme_classic()

cowplot::plot_grid(freq, endb.anno@ggplot + xlim(mut.7.142465001.dt[1, c(2,3)]), ncol=1, align = "v")

### Merged neighboring windows
window.1 <- hotspots.list.merged.ordered[1] %>% makeGRangesFromDataFrame(., ignore.strand = TRUE,
                                                                         start.field="POS",
                                                                         end.field="END_POS",
                                                                         strand.field="CHROM",
                                                                         keep.extra.columns = TRUE)
window.1.dt <- as.data.table(window.1)
freq <- ggplot(data = window.1.dt[c(2:.N)], aes(x=start, fill = region_type)) +
  geom_bar(width = 25, position = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"),
                    name = "Variant type",
                    labels = c("exonic", "intergenic", "intronic", "splice site")) +
  ylab("Variant Frequency") +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(window.1.dt[1, start]-1, window.1.dt[1, end]+1),
                  # here, max(table()) gets the computed frequency value by geom_bar need for upper y_lim
                  ylim = c(0, max(table(window.1.dt[, start]))), 
                  clip = "off")

endb.anno <- autoplot(ensdb, 
                      which = window.1, 
                      columns = c("tx_biotype", "symbol"),
                      names.expr = "symbol:tx_biotype")@ggplot + 
  theme_classic() +
  xlim(window.1.dt[1, start]-1, window.1.dt[1, end]) +
  scale_x_continuous(n.breaks = length(unique(window.1.dt[, start])), 
                     breaks = seq(nearest_lower_limit(window.1.dt[1, start]-1, 5000), 
                                  nearest_upper_limit(window.1.dt[1, end], 5000), 
                                  by = 5000)) +
  xlab(NULL) +
  ylab("Ensembl Transcripts")

cowplot::plot_grid(freq, endb.anno, ncol=1, align = "v", axis = "b")

## TODO
## - Align the axes from top and bottom plots when the legend is on
## - Loop to plot the other hotspots
