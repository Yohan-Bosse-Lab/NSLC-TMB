
#install.packages('VennDiagram')
library(VennDiagram)
# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

library(glue)
library(data.table)

# Creating Venn diagrams from the unfiltered data (before removing variants not passing filters set and using variants not annotated)

# RStudio execution
#output_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#message(paste("Output directory :", output_dir))

# Command line execution (Rscript). Example : Rscript 2-venn_diag_unfiltered.R 0001. Change patient ID 0001 for the desired patient or make a loop to compute all patients.
patient = commandArgs(trailingOnly=TRUE)[1]

if(!is.na(patient)) {
  # If the script is executed from command line, include the following line:
  setwd(glue("/mnt/raid5_hdd/ruelou01/NSLC_{patient}"))

  
  SNPs_mutect <- fread(glue("NSLC_{patient}/SNPs/TNsnv_SNP.txt"), header = FALSE, sep = ":", select = c(2,3,4,5), stringsAsFactors=FALSE)
  SNPs_mutect2 <- fread(glue("NSLC_{patient}/SNPs/Mutect2_SNP.txt"), header = FALSE, sep = ":", select = c(2,3,4,5), stringsAsFactors=FALSE)
  SNPs_TNscope <- fread(glue("NSLC_{patient}/SNPs/TNscope_SNP.txt"), header = FALSE, sep = ":", select = c(2,3,4,5), stringsAsFactors=FALSE)
  
  indels_mutect2 <- fread(glue("NSLC_{patient}/Indels/Mutect2_Ind.txt"), header = FALSE, sep = ":", select = c(2,3,4,5), stringsAsFactors=FALSE)
  indels_TNscope <- fread(glue("NSLC_{patient}/Indels/TNscope_Ind.txt"), header = FALSE, sep = ":", select = c(2,3,4,5), stringsAsFactors=FALSE)
  indels_strelka <- fread(glue("NSLC_{patient}/Indels/strelka_Ind.txt"), header = FALSE, sep = ":", select = c(2,3,4,5), stringsAsFactors=FALSE)
  
  venn.diagram(x=list(SNPs_mutect, SNPs_mutect2, SNPs_TNscope),
               category.names = c("Mutect", "Mutect2", "TNscope"),
               filename = glue('NSLC_{patient}/NSLC_{patient}-unfiltered_Venn_diag_SNP.png'),
               main = "SNPs")
  
  venn.diagram(x=list(indels_mutect2, indels_TNscope, indels_strelka),
               category.names = c("Mutect2", "TNscope", "Strelka"),
               filename = glue('NSLC_{patient}/NSLC_{patient}-unfiltered_Venn_diag_indels.png'),
               main = "Indels")
} else {
  message("Error: Missing patient identifier.")
  quit(save="no",status=1,FALSE)
}

