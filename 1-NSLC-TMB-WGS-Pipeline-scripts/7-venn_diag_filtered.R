# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

#install.packages('VennDiagram')
library(venn)
library(glue)
library(data.table)
library(dplyr)
library(stringr)
library(magrittr)

# Creating Venn diagrams from the filtered data (after removing variants not passing filters set and using annotated variants)

# RStudio execution
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

patients_list <- read.csv("../../93_patients_list.csv", colClasses = "character")

# Command line execution (Rscript). Example : Rscript 7-venn_diag_filtered.R 0001. Change patient ID 0001 for the desired patient or make a loop to compute all patients.
#patient = commandArgs(trailingOnly=TRUE)[1]

col.names = c("Patient_ID", 
              "SNPs_TNscope", "SNPs_Mutect2", "SNPs_Mutect2_TNscope", "SNPs_TNsnv", "SNPs_TNsnv_TNscope", "SNPs_TNsnv_Mutect2", "SNPs_all",
              "Ind_Strelka", "Ind_TNscope", "Ind_TNscope_Strelka", "Ind_Mutect2", "Ind_Mutect2_Strelka", "Ind_Mutect2_TNscope", "Ind_all")

df <- read.table(text = "",
                 col.names = col.names)

for (p in patients_list$x) {
  if(!is.na(p)) {
    patient <- p
    # If the script is executed from command line, include the following line:
    #setwd(glue("/mnt/raid5_hdd/ruelou01/NSLC_{patient}"))
    #
    
    SNPs_TNsnv <- fread(glue("../variants_all_patients/{patient}_TNsnv_completeFilter_SNP.txt"), header = FALSE, stringsAsFactors=FALSE)
    colnames(SNPs_TNsnv) <- "TNsnv"
    SNPs_TNsnv <- SNPs_TNsnv %>% filter(str_detect(TNsnv, 'GL00', negate = TRUE))
    
    SNPs_mutect2 <- fread(glue("../variants_all_patients/{patient}_Mutect2_completeFilter_SNP.txt"), header = FALSE, stringsAsFactors=FALSE)
    colnames(SNPs_mutect2) <- "Mutect2"
    SNPs_mutect2 <- SNPs_mutect2 %>% filter(str_detect(Mutect2, 'GL00', negate = TRUE))
    
    SNPs_TNscope <- fread(glue("../variants_all_patients/{patient}_TNscope_completeFilter_SNP.txt"), header = FALSE, stringsAsFactors=FALSE)
    colnames(SNPs_TNscope) <- "TNscope"
    SNPs_TNscope <- SNPs_TNscope %>% filter(str_detect(TNscope, 'GL00', negate = TRUE))
    
    SNVsVenn <- venn(c(SNPs_TNsnv, SNPs_mutect2, SNPs_TNscope))
    SNVs <- SNVsVenn$counts[2:8]
    
    
    Ind_mutect2 <- fread(glue("../variants_all_patients/{patient}_Mutect2_completeFilter_Ind.txt"), header = FALSE, stringsAsFactors=FALSE)
    colnames(Ind_mutect2) <- "Mutect2"
    Ind_mutect2 <- Ind_mutect2 %>% filter(str_detect(Mutect2, 'GL00', negate = TRUE))
    
    Ind_TNscope <- fread(glue("../variants_all_patients/{patient}_TNscope_completeFilter_Ind.txt"), header = FALSE, stringsAsFactors=FALSE)
    colnames(Ind_TNscope) <- "TNscope"
    Ind_TNscope <- Ind_TNscope %>% filter(str_detect(TNscope, 'GL00', negate = TRUE))
    
    Ind_Strelka <- fread(glue("../variants_all_patients/{patient}_Strelka_completeFilter_Ind.txt"), header = FALSE, stringsAsFactors=FALSE)
    colnames(Ind_Strelka) <- "Strelka"
    Ind_Strelka <- Ind_Strelka %>% filter(str_detect(Strelka, 'GL00', negate = TRUE))
    
    
    IndVenn <- venn(c(Ind_mutect2, Ind_TNscope, Ind_Strelka))
    Ind <- IndVenn$counts[2:8]
    
    df[nrow(df) + 1,] <- c(patient, SNVs, Ind)
    
  } else {
    message("Error: Missing patient identifier argument.")
    quit(save="default",status=1,FALSE)
  }
}

df[, c(2:15)] <- sapply(df[, c(2:15)], as.numeric)
colMeans(df[2:15])

apply(df[2:15], 2, median)
apply(df[2:15], 2, min)
apply(df[2:15], 2, max)


