# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

library(glue)
library(data.table)

# Reading vcf files from variant calling tools and applying filters (read depth for normal and tumor samples, allele frequency)

# RStudio execution
#patient = "0001"
#output_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#message(paste("Output directory:", output_dir))

# Command line execution (Rscript)
patient = commandArgs(trailingOnly=TRUE)[1]


if(!is.na(patient)) {
  setwd(glue("/mnt/raid5_hdd/ruelou01/NSLC_{patient}"))
  message(glue("Output directory: /mnt/raid5_hdd/ruelou01/NSLC_{patient}")) 

  # Mutect2
  
  message("Mutect2 data import.")
  Mutect2 <- readLines(glue("Mutect2/{patient}_Mutect2_norm_ID.vcf.gz"))
  Mutect2_data <- read.table(glue("Mutect2/{patient}_Mutect2_norm_ID.vcf.gz"), sep = "\t")
  Mutect2<-Mutect2[-(grep("#CHROM",Mutect2)+1):-(length(Mutect2))]
  Mutect2_names<-unlist(strsplit(Mutect2[length(Mutect2)],"\t"))
  names(Mutect2_data)<-Mutect2_names
  
  message("Mutect2 table creation.")
  for (i in 1L:as.integer(nrow(Mutect2_data))) {
    if(i==1L)
    {
      Mutect2_table <- Mutect2_data[,1:5]
      Mutect2_table["depth_tumor"] <- 0
      Mutect2_table["depth_normal"] <- 0
      Mutect2_table["variant_readCount_tumor"] <-0
      Mutect2_table["variantAF_normal"] <- 0
    }
    tumor_numbers <- as.numeric(unlist(strsplit(as.character(Mutect2_data[i,]$Mutect2_TUMOR), "\\:|\\,"))[c(2,3)])
    normal_numbers <- as.numeric(unlist(strsplit(as.character(Mutect2_data[i,]$Mutect2_NORMAL), "\\:|\\,"))[c(2,3)])
    set(Mutect2_table, i, 6L, sum(tumor_numbers))
    set(Mutect2_table, i, 7L, sum(normal_numbers))
    set(Mutect2_table, i, 8L, tumor_numbers[2])
    set(Mutect2_table, i, 9L, normal_numbers[2]/sum(normal_numbers))
  }
  # Mutect2 filtering
  Mutect2_table <- na.omit(Mutect2_table[Mutect2_table$depth_tumor > 12 & Mutect2_table$depth_normal > 6 & Mutect2_table$variant_readCount_tumor > 5 & Mutect2_table$variantAF_normal < 0.02,])
  
  message("Mutect2 table and filtered vcf export.")
  fwrite(Mutect2_table, file=glue("{patient}_Mutect2.csv"), sep=";", dec = ".", quote = TRUE)
  write(Mutect2, file=glue("Mutect2/{patient}_Mutect2_preFilter.vcf"))
  write.table(merge(Mutect2_data, Mutect2_table, by="ID", sort = FALSE)[,c(2,3,1,4,5:11)], file=glue("Mutect2/{patient}_Mutect2_preFilter.vcf"), append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  
  # TNscope
  
  message("TNscope data import.")
  TNscope <- readLines(glue("TNscope/{patient}_TNscope_norm_ID.vcf.gz"))
  TNscope_data <- read.table(glue("TNscope/{patient}_TNscope_norm_ID.vcf.gz"), sep = "\t")
  TNscope<-TNscope[-(grep("#CHROM",TNscope)+1):-(length(TNscope))]
  TNscope_names<-unlist(strsplit(TNscope[length(TNscope)],"\t"))
  names(TNscope_data)<-TNscope_names
  
  message("TNscope table creation.")
  TNscope_range <- as.integer(rownames(TNscope_data[TNscope_data$FILTER == "PASS",]))
  index=1L
  for (i in TNscope_range) {
    if(index==1L){
      TNscope_table <- TNscope_data[TNscope_data$FILTER == "PASS",1:5]
      TNscope_table["depth_tumor"] <- 0
      TNscope_table["depth_normal"] <- 0
      TNscope_table["variant_readCount_tumor"] <-0
      TNscope_table["variantAF_normal"] <- 0
    }
    tumor_numbers <- suppressWarnings(as.numeric(unlist(strsplit(as.character(TNscope_data[i,]$TNscope_TUMOR), "\\:|\\,"))[c(2,3)]))
    normal_numbers <- suppressWarnings(as.numeric(unlist(strsplit(as.character(TNscope_data[i,]$TNscope_NORMAL), "\\:|\\,"))[c(2,3)]))
    set(TNscope_table, index, 6L, sum(tumor_numbers))
    set(TNscope_table, index, 7L, sum(normal_numbers))
    set(TNscope_table, index, 8L, tumor_numbers[2])
    set(TNscope_table, index, 9L, normal_numbers[2]/sum(normal_numbers))
    index=index+1L
  }
  # TNscope filtering
  TNscope_table <- na.omit(TNscope_table[TNscope_table$depth_tumor > 12 & TNscope_table$depth_normal > 6 & TNscope_table$variant_readCount_tumor > 5 & TNscope_table$variantAF_normal < 0.02,])
  
  message("TNscope table and filtered vcf export.")
  fwrite(TNscope_table, file=glue("{patient}_TNscope.csv"), sep=";", dec = ".", quote = TRUE)
  write(TNscope, file=glue("TNscope/{patient}_TNscope_preFilter.vcf"))
  write.table(merge(TNscope_data, TNscope_table, by="ID", sort=FALSE)[,c(2,3,1,4:11)], file=glue("TNscope/{patient}_TNscope_preFilter.vcf"), append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  
  # Strelka
  
  message("Strelka data import.")
  Strelka <- readLines(glue("Strelka/{patient}_Strelka_norm_ID.indels.vcf.gz"))
  Strelka_data <- read.table(glue("Strelka/{patient}_Strelka_norm_ID.indels.vcf.gz"), sep = "\t")
  Strelka<-Strelka[-(grep("#CHROM",Strelka)+1):-(length(Strelka))]
  Strelka_names<-unlist(strsplit(Strelka[length(Strelka)],"\t"))
  names(Strelka_data)<-Strelka_names
  
  message("Strelka table creation.")
  Strelka_range <- as.integer(rownames(Strelka_data[Strelka_data$FILTER == "PASS",]))
  index=1L
  for (i in Strelka_range) {
    if(index==1L)
    {
      Strelka_table <- Strelka_data[Strelka_data$FILTER == "PASS",1:5]
      Strelka_table["depth_tumor"] <- 0
      Strelka_table["depth_normal"] <- 0
      Strelka_table["variant_readCount_tumor"] <-0
      Strelka_table["variantAF_normal"] <- 0
    }
    # Strelka output format : first comma-delimited value from TAR (3rd column) is the REF counts, first comma-delimited value from TIR (5th column) is the ALT counts.
    tumor_numbers <- as.numeric(unlist(strsplit(as.character(Strelka_data[i,]$Strelka_TUMOR), "\\:|\\,"))[c(3,5)])
    normal_numbers <- as.numeric(unlist(strsplit(as.character(Strelka_data[i,]$Strelka_NORMAL), "\\:|\\,"))[c(3,5)])
    set(Strelka_table, index, 6L, sum(tumor_numbers))
    set(Strelka_table, index, 7L, sum(normal_numbers))
    set(Strelka_table, index, 8L, tumor_numbers[2])
    set(Strelka_table, index, 9L, normal_numbers[2]/sum(normal_numbers))
    index=index+1L
  }
  # Strelka filtering
  Strelka_table <- na.omit(Strelka_table[Strelka_table$depth_tumor > 12 & Strelka_table$depth_normal > 6 & Strelka_table$variant_readCount_tumor > 5 & Strelka_table$variantAF_normal < 0.02,])
  
  message("Strelka table and filtered vcf export.")
  fwrite(Strelka_table, file=glue("{patient}_Strelka.csv"), sep=";", dec = ".", quote = TRUE)
  write(Strelka, file=glue("Strelka/{patient}_Strelka_preFilter.vcf"))
  write.table(merge(Strelka_data, Strelka_table, by="ID", sort=FALSE)[,c(2,3,1,4:11)], file=glue("Strelka/{patient}_Strelka_preFilter.vcf"), append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  
  # Mutect (TNsnv)
  
  message("Mutect data import.")
  Mutect <- readLines(glue("TNsnv/{patient}_TNsnv_norm_ID.vcf.gz"))
  Mutect_data <- read.table(glue("TNsnv/{patient}_TNsnv_norm_ID.vcf.gz"), sep = "\t")
  Mutect<-Mutect[-(grep("#CHROM",Mutect)+1):-(length(Mutect))]
  Mutect_names<-unlist(strsplit(Mutect[length(Mutect)],"\t"))
  names(Mutect_data)<-Mutect_names
  
  message("Mutect table creation.")
  Mutect_range <- as.integer(rownames(Mutect_data[Mutect_data$FILTER == "PASS",]))
  index=1L
  for (i in Mutect_range) {
    if(index==1L){
      Mutect_table <- Mutect_data[Mutect_data$FILTER == "PASS",1:5]
      Mutect_table["depth_tumor"] <- 0
      Mutect_table["depth_normal"] <- 0
      Mutect_table["variant_readCount_tumor"] <-0
      Mutect_table["variantAF_normal"] <- 0
    }
    tumor_numbers <- as.numeric(unlist(strsplit(as.character(Mutect_data[i,]$TNsnv_TUMOR), "\\:|\\,"))[c(2,3)])
    normal_numbers <- as.numeric(unlist(strsplit(as.character(Mutect_data[i,]$TNsnv_NORMAL), "\\:|\\,"))[c(2,3)])
    set(Mutect_table, index, 6L, sum(tumor_numbers))
    set(Mutect_table, index, 7L, sum(normal_numbers))
    set(Mutect_table, index, 8L, tumor_numbers[2])
    set(Mutect_table, index, 9L, normal_numbers[2]/sum(normal_numbers))
    index=index+1L
  }
  # Mutect filtering
  Mutect_table <- na.omit(Mutect_table[Mutect_table$depth_tumor > 12 & Mutect_table$depth_normal > 6 & Mutect_table$variant_readCount_tumor > 5 & Mutect_table$variantAF_normal < 0.02,])
  
  message("Mutect table and filtered vcf export.")
  fwrite(Mutect_table, file=glue("{patient}_TNsnv.csv"), sep=";", dec = ".", quote = TRUE)
  write(Mutect, file=glue("TNsnv/{patient}_TNsnv_preFilter.vcf"))
  write.table(merge(Mutect_data, Mutect_table, by="ID", sort=FALSE)[,c(2,3,1,4:11)], file=glue("TNsnv/{patient}_TNsnv_preFilter.vcf"), append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

} else {
  message("Error: Missing patient identifier argument.")
  quit(save="default",status=1,FALSE)
}



