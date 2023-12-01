
## TODO
## - Compare TMBs with SBS signatures
## - Make a list of cancer protectors, cancer drivers and DNA damage repair genes to find mutations

library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(biomaRt)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/n82_NSLC_surv_data.xlsx", sheetIndex = 1))

# Fetching and merging all patients variants
VEP_data_all_patients <- fread("Data/n82_NSLC_variants_data.csv", keepLeadingZeros = TRUE)

# DNA damage response genes
DDR.genes <- setDT(readxl::read_excel("Data/NIHMS962902_DDR_genes.xlsx", range= "A4:AC280"))
DDR.gene.mutations <- setNames(lapply(DDR.genes[,`Gene Symbol`], function(gene) VEP_data_all_patients[SYMBOL == gene]), 
                               nm= DDR.genes[,`Gene Symbol`])


lapply(names(DDR.gene.mutations), function(x) mutate(surv.dt, ))

################################################################################
#################################### Archive ###################################
################################################################################

# Searching for mutated DDR genes by genome coordinates
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
DDR.genes.info  <- setDT(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                               filters=c('entrezgene_id', 'hgnc_symbol'),
                               values=list(DDR.genes[, `Entrez Gene ID`], DDR.genes[, `Gene Symbol`]),
                               mart=ensembl))
# Order chromosome names
DDR.genes.info <- DDR.genes.info[str_order(chromosome_name, numeric = TRUE),]
# Remove non standard chromosome names
DDR.genes.info <- DDR.genes.info[chromosome_name %in% c(1:22, 'X', 'Y')]

DDR.genes.mutated <-  setNames(lapply(1:nrow(DDR.genes.info), function(x) VEP_data_all_patients[POS >= DDR.genes.info[x, start_position] &
                                                                                                  END_POS <= DDR.genes.info[x, end_position] &
                                                                                                  CHROM == DDR.genes.info[x, chromosome_name]]),
                               nm = DDR.genes.info[, hgnc_symbol])


