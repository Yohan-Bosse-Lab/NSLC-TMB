if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("org.Hs.eg.db")

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Another attempt at exome length, but not not the right result.
exome <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
exome_ranges <- IRangesList(exome)
exome_ranges$chr1@start
sum(unlist(width(exome_ranges)[1:24]))

cds <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
cds_ranges <- IRangesList(cds)
sum(unlist(width(cds_ranges)[1:24]))

# Gene length : from https://www.biostars.org/p/62583/#62634
hg19GeneLengths <- function(symbols)
{
  require(TxDb.Hsapiens.UCSC.hg19.knownGene) 
  require(org.Hs.eg.db)
  cds.db = cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')    
  egs = unlist(mget(symbols[symbols %in% keys(org.Hs.egSYMBOL2EG)], org.Hs.egSYMBOL2EG))
  sapply(egs,function(eg)
  {
    cds = cds.db[[eg]]
    cds = reduce(cds)
    sum(width(cds))
  })
}
F1CDx <- hg19GeneLengths(c("ABL1","ACVR1B","AKT1","AKT2","AKT3","ALK","ALOX12B","AMER1","APC","AR","ARAF","ARFRP1","ARID1A","ASXL1","ATM","ATR","ATRX","AURKA","AURKB","AXIN1","AXL","BAP1","BARD1","BCL2","BCL2L1","BCL2L2","BCL6","BCOR","BCORL1","BRAF","BRCA1","BRCA2","BRD4","BRIP1","BTG1","BTG2","BTK","C11orf30","CALR","CARD11","CASP8","CBFB","CBL","CCND1","CCND2","CCND3","CCNE1","CD22","CD274","CD70","CD79A","CD79B","CDC73","CDH1","CDK12","CDK4","CDK6","CDK8","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CEBPA","CHEK1","CHEK2","CIC","CREBBP","CRKL","CSF1R","CSF3R","CTCF","CTNNA1","CTNNB1","CUL3","CUL4A","CXCR4","CYP17A1","DAXX","DDR1","DDR2","DIS3","DNMT3A","DOT1L","EED","EGFR","EP300","EPHA3","EPHB1","EPHB4","ERBB2","ERBB3","ERBB4","ERCC4","ERG","ERRFI1","ESR1","EZH2","FAM46C","FANCA","FANCC","FANCG","FANCL","FAS","FBXW7","FGF10","FGF12","FGF14","FGF19","FGF23","FGF3","FGF4","FGF6","FGFR1","FGFR2","FGFR3","FGFR4","FH","FLCN","FLT1","FLT3","FOXL2","FUBP1","GABRA6","GATA3","GATA4","GATA6","GID4","GNA11","GNA13","GNAQ","GNAS","GRM3","GSK3B","H3F3A","HDAC1","HGF","HNF1A","HRAS","HSD3B1","ID3","IDH1","IDH2","IGF1R","IKBKE","IKZF1","INPP4B","IRF2","IRF4","IRS2","JAK1","JAK2","JAK3","JUN","KDM5A","KDM5C","KDM6A","KDR","KEAP1","KEL","KIT","KLHL6","KMT2A","KMT2D","KRAS","LTK","LYN","MAF","MAP2K1","MAP2K2","MAP2K4","MAP3K1","MAP3K13","MAPK1","MCL1","MDM2","MDM4","MED12","MEF2B","MEN1","MERTK","MET","MITF","MKNK1","MLH1","MPL","MRE11A","MSH2","MSH3","MSH6","MST1R","MTAP","MTOR","MUTYH","MYC","MYCL","MYCN","MYD88","NBN","NF1","NF2","NFE2L2","NFKBIA","NKX2-1","NOTCH1","NOTCH2","NOTCH3","NPM1","NRAS","NT5C2","NTRK1","NTRK2","NTRK3","P2RY8","PALB2","PARK2","PARP1","PARP2","PARP3","PAX5","PBRM1","PDCD1","PDCD1LG2","PDGFRA","PDGFRB","PDK1","PIK3C2B","PIK3C2G","PIK3CA","PIK3CB","PIK3R1","PIM1","PMS2","POLD1","POLE","PPARG","PPP2R1A","PPP2R2A","PRDM1","PRKAR1A","PRKCI","PTCH1","PTEN","PTPN11","PTPRO","QKI","RAC1","RAD21","RAD51","RAD51B","RAD51C","RAD51D","RAD52","RAD54L","RAF1","RARA","RB1","RBM10","REL","RET","RICTOR","RNF43","ROS1","RPTOR","SDHA","SDHB","SDHC","SDHD","SETD2","SF3B1","SGK1","SMAD2","SMAD4","SMARCA4","SMARCB1","SMO","SNCAIP","SOCS1","SOX2","SOX9","SPEN","SPOP","SRC","STAG2","STAT3","STK11","SUFU","SYK","TBX3","TEK","TET2","TGFBR2","TIPARP","TNFAIP3","TNFRSF14","TP53","TSC1","TSC2","TYRO3","U2AF1","VEGFA","VHL","WHSC1","WHSC1L1","WT1","XPO1","XRCC2","ZNF217","ZNF703"))
sum(F1CDx)

Illu500 <- hg19GeneLengths(c("ABL1","ABL2","ACVR1","ACVR1B","AKT1","AKT2","AKT3","ALK","ALOX12B","ANKRD11","ANKRD26","APC","AR","ARAF","ARFRP1","ARID1A","ARID1B","ARID2","ARID5B","ASXL1","ASXL2","ATM","ATR","ATRX","AURKA","AURKB","AXIN1","AXIN2","AXL","B2M","BAP1","BARD1","BBC3","BCL10","BCL2","BCL2L1","BCL2L11","BCL2L2","BCL6","BCOR","BCORL1","BCR","BIRC3","BLM","BMPR1A","BRAF","BRCA1","BRCA2","BRD4","BRIP1","BTG1","BTK","C11orf30","CALR","CARD11","CASP8","CBFB","CBL","CCND1","CCND2","CCND3","CCNE1","CD274","CD276","CD74","CD79A","CD79B","CDC73","CDH1","CDK12","CDK4","CDK6","CDK8","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CEBPA","CENPA","CHD2","CHD4","CHEK1","CHEK2","CIC","CREBBP","CRKL","CRLF2","CSF1R","CSF3R","CSNK1A1","CTCF","CTLA4","CTNNA1","CTNNB1","CUL3","CUX1","CXCR4","CYLD","DAXX","DCUN1D1","DDR2","DDX41","DHX15","DICER1","DIS3","DNAJB1","DNMT1","DNMT3A","DNMT3B","DOT1L","E2F3","EED","EGFL7","EGFR","EIF1AX","EIF4A2","EIF4E","EML4","EP300","EPCAM","EPHA3","EPHA5","EPHA7","EPHB1","ERBB2","ERBB3","ERBB4","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERG","ERRFI1","ESR1","ETS1","ETV1","ETV4","ETV5","ETV6","EWSR1","EZH2","FAM123B","FAM175A","FAM46C","FANCA","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI","FANCL","FAS","FAT1","FBXW7","FGF1","FGF10","FGF14","FGF19","FGF2","FGF23","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGFR1","FGFR2","FGFR3","FGFR4","FH","FLCN","FLI1","FLT1","FLT3","FLT4","FOXA1","FOXL2","FOXO1","FOXP1","FRS2","FUBP1","FYN","GABRA6","GATA1","GATA2","GATA3","GATA4","GATA6","GEN1","GID4","GLI1","GNA11","GNA13","GNAQ","GNAS","GPR124","GPS2","GREM1","GRIN2A","GRM3","GSK3B","H3F3A","H3F3B","H3F3C","HGF","HIST1H1C","HIST1H2BD","HIST1H3A","HIST1H3B","HIST1H3C","HIST1H3D","HIST1H3E","HIST1H3F","HIST1H3G","HIST1H3H","HIST1H3I","HIST1H3J","HIST2H3A","HIST2H3C","HIST2H3D","HIST3H3","HLA-A","HLA-B","HLA-C","HNF1A","HNRNPK","HOXB13","HRAS","HSD3B1","HSP90AA1","ICOSLG","ID3","IDH1","IDH2","IFNGR1","IGF1","IGF1R","IGF2","IKBKE","IKZF1","IL10","IL7R","INHA","INHBA","INPP4A","INPP4B","INSR","IRF2","IRF4","IRS1","IRS2","JAK1","JAK2","JAK3","JUN","KAT6A","KDM5A","KDM5C","KDM6A","KDR","KEAP1","KEL","KIF5B","KIT","KLF4","KLHL6","KMT2B","KMT2C","KMT2D","KRAS","LAMP1","LATS1","LATS2","LMO1","LRP1B","LYN","LZTR1","MAGI2","MALT1","MAP2K1","MAP2K2","MAP2K4","MAP3K1","MAP3K13","MAP3K14","MAP3K4","MAPK1","MAPK3","MAX","MCL1","MDC1","MDM2","MDM4","MED12","MEF2B","MEN1",
                             "MET","MGA","MITF","MLH1","MLL","MLLT3","MPL","MRE11A","MSH2","MSH3","MSH6","MST1","MST1R","MTOR","MUTYH","MYB","MYC","MYCL1","MYCN","MYD88","MYOD1","NAB2","NBN","NCOA3","NCOR1","NEGR1","NF1","NF2","NFE2L2","NFKBIA","NKX2-1","NKX3-1","NOTCH1","NOTCH2","NOTCH3","NOTCH4","NPM1","NRAS","NRG1","NSD1","NTRK1","NTRK2","NTRK3","NUP93","NUTM1","PAK1","PAK3","PAK7","PALB2","PARK2","PARP1","PAX3","PAX5","PAX7","PAX8","PBRM1","PDCD1","PDCD1LG2","PDGFRA","PDGFRB","PDK1","PDPK1","PGR","PHF6","PHOX2B","PIK3C2B","PIK3C2G","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIM1","PLCG2","PLK2","PMAIP1","PMS1","PMS2","PNRC1","POLD1","POLE","PPARG","PPM1D","PPP2R1A","PPP2R2A","PPP6C","PRDM1","PREX2","PRKAR1A","PRKCI","PRKDC","PRSS8","PTCH1","PTEN","PTPN11","PTPRD","PTPRS","PTPRT","QKI","RAB35","RAC1","RAD21","RAD50","RAD51","RAD51B","RAD51C","RAD51D","RAD52","RAD54L","RAF1","RANBP2","RARA","RASA1","RB1","RBM10","RECQL4","REL","RET","RFWD2","RHEB","RHOA","RICTOR","RIT1","RNF43","ROS1","RPS6KA4","RPS6KB1","RPS6KB2","RPTOR","RUNX1","RUNX1T1","RYBP","SDHA","SDHAF2","SDHB","SDHC","SDHD","SETBP1","SETD2","SF3B1","SH2B3","SH2D1A","SHQ1","SLIT2","SLX4","SMAD2","SMAD3","SMAD4","SMARCA4","SMARCB1","SMARCD1","SMC1A","SMC3","SMO","SNCAIP","SOCS1","SOX10","SOX17","SOX2","SOX9","SPEN","SPOP","SPTA1","SRC","SRSF2","STAG1","STAG2","STAT3","STAT4","STAT5A","STAT5B","STK11","STK40","SUFU","SUZ12","SYK","TAF1","TBX3","TCEB1","TCF3","TCF7L2","TERC","TERT","TET1","TET2","TFE3","TFRC","TGFBR1","TGFBR2","TMEM127","TMPRSS2","TNFAIP3","TNFRSF14","TOP1","TOP2A","TP53","TP63","TRAF2","TRAF7","TSC1","TSC2","TSHR","U2AF1","VEGFA","VHL","VTCN1","WISP3","WT1","XIAP","XPO1","XRCC2","YAP1","YES1","ZBTB2","ZBTB7A","ZFHX3","ZNF217","ZNF703","ZRSR2"))
sum(Illu500)

MSKI468 <- hg19GeneLengths(c("ABL1","ACVR1","AGO2","AKT1","AKT2","AKT3","ALK","ALOX12B","AMER1","ANKRD11","APC","AR","ARAF","ARID1A","ARID1B","ARID2","ARID5B","ASXL1","ASXL2","ATM","ATR","ATRX","AURKA","AURKB","AXIN1","AXIN2","AXL","B2M","BABAM1","BAP1","BARD1","BBC3","BCL10","BCL2","BCL2L1","BCL2L11","BCL6","BCOR","BIRC3","BLM","BMPR1A","BRAF","BRCA1","BRCA2","BRD4","BRIP1","BTK","CALR","CARD11","CARM1","CASP8","CBFB","CBL","CCND1","CCND2","CCND3","CCNE1","CD274","CD276","CD79A","CD79B","CDC42","CDC73","CDH1","CDK12","CDK4","CDK6","CDK8","CDKN1A","CDKN1B","CDKN2Ap14ARF","CDKN2Ap16INK4A","CDKN2B","CDKN2C","CEBPA","CENPA","CHEK1","CHEK2","CIC","CREBBP","CRKL","CRLF2","CSDE1","CSF1R","CSF3R","CTCF","CTLA4","CTNNB1","CUL3","CXCR4","CYLD","CYSLTR2","DAXX","DCUN1D1","DDR2","DICER1","DIS3","DNAJB1","DNMT1","DNMT3A","DNMT3B","DOT1L","DROSHA","DUSP4","E2F3","EED","EGFL7","EGFR","EIF1AX","EIF4A2","EIF4E","ELF3","EP300","EPAS1","EPCAM","EPHA3","EPHA5","EPHA7","EPHB1","ERBB2","ERBB3","ERBB4","ERCC2","ERCC3","ERCC4","ERCC5","ERF","ERG","ERRFI1","ESR1","ETV1","ETV6","EZH1","EZH2","FAM175A","FAM46C","FAM58A","FANCA","FANCC","FAT1","FBXW7","FGF19","FGF3","FGF4","FGFR1","FGFR2","FGFR3","FGFR4","FH","FLCN","FLT1","FLT3","FLT4","FOXA1","FOXL2","FOXO1","FOXP1","FUBP1","FYN","GATA1","GATA2","GATA3","GLI1","GNA11","GNAQ","GNAS","GPS2","GREM1","GRIN2A","GSK3B","H3F3A","H3F3B","H3F3C","HGF","HIST1H1C","HIST1H2BD","HIST1H3A","HIST1H3B","HIST1H3C","HIST1H3D","HIST1H3E","HIST1H3F","HIST1H3G","HIST1H3H","HIST1H3I","HIST1H3J","HIST2H3C","HIST2H3D","HIST3H3","HLA-A","HLA-B","HNF1A","HOXB13","HRAS","ICOSLG","ID3","IDH1","IDH2","IFNGR1","IGF1","IGF1R","IGF2","IKBKE","IKZF1","IL10","IL7R","INHA","INHBA","INPP4A","INPP4B","INPPL1","INSR","IRF4","IRS1","IRS2","JAK1","JAK2","JAK3","JUN","KDM5A","KDM5C","KDM6A","KDR","KEAP1","KIT","KLF4","KMT2A","KMT2B","KMT2C","KMT2D","KNSTRN","KRAS","LATS1","LATS2","LMO1","LYN","MALT1","MAP2K1","MAP2K2","MAP2K4","MAP3K1","MAP3K13","MAP3K14","MAPK1","MAPK3","MAPKAP1","MAX","MCL1","MDC1","MDM2","MDM4","MED12","MEF2B","MEN1","MET","MGA","MITF","MLH1","MPL","MRE11A","MSH2","MSH3","MSH6","MSI1","MSI2","MST1","MST1R","MTOR","MUTYH","MYC","MYCL1","MYCN","MYD88","MYOD1","NBN","NCOA3","NCOR1","NEGR1","NF1","NF2","NFE2L2","NFKBIA","NKX2-1","NKX3-1","NOTCH1","NOTCH2","NOTCH3","NOTCH4","NPM1","NRAS","NSD1","NTHL1","NTRK1","NTRK2","NTRK3","NUF2","NUP93","PAK1","PAK7","PALB2","PARK2","PARP1","PAX5","PBRM1","PDCD1","PDCD1LG2","PDGFRA","PDGFRB","PDPK1","PGR","PHOX2B","PIK3C2G","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIM1","PLCG2","PLK2","PMAIP1","PMS1","PMS2","PNRC1","POLD1","POLE","PPARG","PPM1D","PPP2R1A","PPP4R2","PPP6C","PRDM1","PRDM14","PREX2","PRKAR1A","PRKCI","PRKD1","PTCH1","PTEN","PTP4A1","PTPN11","PTPRD","PTPRS","PTPRT","RAB35","RAC1","RAC2","RAD21","RAD50","RAD51","RAD51B","RAD51C","RAD51D","RAD52","RAD54L","RAF1","RARA","RASA1","RB1","RBM10","RECQL","RECQL4","REL","RET","RFWD2","RHEB","RHOA","RICTOR","RIT1","RNF43","ROS1","RPS6KA4","RPS6KB2","RPTOR","RRAGC","RRAS","RRAS2","RTEL1","RUNX1","RXRA","RYBP","SDHA","SDHAF2","SDHB","SDHC","SDHD","SESN1","SESN2","SESN3","SETD2","SETD8","SF3B1","SH2B3","SH2D1A","SHOC2","SHQ1","SLX4","SMAD2","SMAD3","SMAD4","SMARCA4","SMARCB1","SMARCD1","SMO","SMYD3","SOCS1","SOS1","SOX17","SOX2","SOX9","SPEN","SPOP","SPRED1","SRC","SRSF2","STAG2","STAT3","STAT5A","STAT5B","STK11","STK19","STK40","SUFU","SUZ12","SYK","TAP1","TAP2","TBX3","TCEB1","TCF3","TCF7L2","TEK","TERT","TET1","TET2","TGFBR1","TGFBR2","TMEM127","TMPRSS2","TNFAIP3","TNFRSF14","TOP1","TP53","TP53BP1","TP63","TRAF2","TRAF7","TSC1","TSC2","TSHR","U2AF1","UPF1","VEGFA","VHL","VTCN1","WHSC1","WHSC1L1","WT1","WWTR1","XIAP","XPO1","XRCC2","YAP1","YES1","ZFHX3"))
sum(MSKI468)

NeoPlus <- hg19GeneLengths(c("ABL1","ABL2","ACVR1B","AKT1","AKT2","AKT3","ALK","ALMS1","AMER1","APC","APLNR","AR","ARAF","ARHGEF12","ARID1A","ARID1B","ARID2","ASXL1","ATAD5","ATM","ATR","ATRX","AURKA","AURKB","AXIN2","AXL","B2M","BAP1","BARD1","BCL2","BCL2L1","BCL6","BLM","BMPR1A","BMS1","BRAF","BRCA1","BRCA2","BRD2","BRD3","BRD4","BRIP1","BTK","BUB1B","CARD11","CASP8","CBFB","CBL","CCND1","CCND2","CCNE1","CD274","CD58","CD74","CD79A","CD79B","CDC73","CDH1","CDK12","CDK4","CDK6","CDK8","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CEBPA","CHD4","CHEK1","CHEK2","CIC","CLSPN","CREBBP","CRKL","CSF1R","CTCF","CTNNA1","CTNNB1","CUL3","DAXX","DCUN1D4","DDB2","DDR2","DICER1","DOT1L","EGFR","EML4","EMSY","EP300","EPHA3","EPHA5","EPHA7","EPHB1","ERBB2","ERBB3","ERBB4","ERCC3","ERCC4","ERCC5","ERCC6","ERRFI1","ESR1","ETV6","EZH2","FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI","FANCL","FANCM","FAS","FAT1","FBXW7","FGF19","FGF23","FGF3","FGF4","FGF6","FGFR1","FGFR2","FGFR3","FGFR4","FLCN","FLT1","FLT3","FOXO3","FOXP1","FRS2","GATA1","GEN1","GLI1","GNA11","GNA13","GNAI2","GNAQ","GNAS","GNAT2","GSK3B","H3F3A","H3F3B","HDAC2","HGF","HLTF","HRAS","HSP90AA1","IDH1","IDH2","IFNGR1","IFNGR2","IGF1R","IGF2","IGF2R","IKBKB","IKBKE","IL21R","INHBA","INPP4B","IRF1","IRF2","IRF4","JAK1","JAK2","JAK3","JUN","KAT6A","KDR","KEAP1","KIF5B","KIT","KMT2A","KMT2B","KMT2C","KMT2D","KRAS","KSR1","LRP1B","LYN","MAD2L1","MAGI2","MAP2K1","MAP2K3","MAP2K4","MAP3K1","MCL1","MDM2","MDM4","MED12","MEN1","MET","MITF","MLANA","MLH1","MLH3","MPL","MRE11","MSH2","MSH3","MSH5","MSH6","MST1R","MTOR","MUTYH","MYC","MYCL","MYCN","NBN","NCOA3","NCOA4","NF1","NF2","NFE2L2","NFKBIA","NOTCH2","NOTCH4","NPM1","NRAS","NRG1","NSD1","NTRK1","NTRK2","NTRK3","PALB2","PBRM1","PDCD1LG2","PDGFRA","PDGFRB","PDK1","PIK3C2B","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PLCG2","PMS2","POLD1","POLE","POLG","POLH","POT1","PPM1D","PPP2R1A","PRDM1","PREX2","PRKCI","PRKDC","PRKN","PSMB5","PTCH1","PTEN","PTPN11","PTPRC","PTPRK","PTPRT","QKI","RAC1","RAD50","RAD51","RAD51B","RAD51C","RAD51D","RAD54L","RAD54L2","RAF1","RANBP2","RARA","RB1","RECQL4","RET","REV3L","RICTOR","RIT1","RNF43","RNPS1","ROS1","RPL10A","RPL23","RPTOR","SDHA","SDHAF2","SDHB","SDHC","SDHD","SERPINB3","SERPINB4","SETD1A","SETD1B","SETD2","SF3B1","SLIT2","SLX4","SMAD2","SMAD3","SMAD4","SMARCA4","SMARCB1","SMO","SOX2","SPEN","SPOP","SRP54","STAG2","STAT3","STK11","SUFU","SYK","TAF3","TAP1","TAP2","TAPBP","TCF7L2","TGFBR2","TNFAIP3","TNFRSF14","TOP1","TOP2A","TP53","TP53BP1","TRRAP","TSC1","TSC2","TSHR","U2AF1","U2AF1L5","VEGFA","VHL","WISP3","WRN","WT1","XPA","XPC","XRCC1","ZFHX3","ZNF217"))
sum(NeoPlus)

OncoMineTMB <- hg19GeneLengths(c("ABL2","ACVR2A","ADAMTS20","AFF1","AFF3","AKAP9","APC","ARID2","ARNT","ATF1","AURKA","AURKB","AURKC","BAI3","BCL10","BCL11A","BCL11B","BCL2","BCL2L1","BCL2L2","BCL3","BCL6","BCL9","BCR","BIRC2","BIRC3","BIRC5","BLM","BLNK","BMPR1A","BRD3","BTK","BUB1B","CARD11","CASC5","CCND2","CCNE1","CD79A","CD79B","CDC73","CDH1","CDH11","CDH2","CDH20","CDH5","CDK8","CDKN2C","CIC","CKS1B","CMPK1","COL1A1","CRBN","CREB1","CRKL","CRTC1","CSMD3","CTNNA1","CTNNB1","CYLD","CYP2C19","CYP2D6","DAXX","DCC","DDB2","DDIT3","DEK","DICER1","DPYD","DST","EML4","EP300","EP400","EPHA3","EPHA7","EPHB1","EPHB4","EPHB6","ERCC1","ERCC3","ERCC4","ERCC5","ERG","ETS1","ETV1","ETV4","EXT1","EXT2","FAM123B","FANCC","FANCF","FANCG","FANCJ","FAS","FH","FLCN","FLl1","FLT1","FLT4","FN1","FOXL2","FOXO1","FOXO3","FOXP1","FOXP4","FZR1","G6PD","GATA1","GATA2","GATA3","GDNF","GPR124","GRM8","GUCY1A2","HCAR1","HIF1A","HLF","HOOK3","HSP90AA1","HSP90AB1","ICK","IGF1R","IGF2","IGF2R","IKBKB","IKBKE","IKZF1","IL2","IL21R","IL6ST","IL7R","ING4","IRF4","IRS2","ITGA10","ITGA9","ITGB2","ITGB3","JAK1","JAK3","JUN","KAT6A","KAT6B","KDM5C","KDM6A","KEAP1","KLF6","LAMP1","LCK","LIFR","LPHN3","LPP","LRP1B","LTF","LTK","MAF","MAFB","MAGEA1","MAGl1","MALT1","MAML2","MAP3K7","MAPK8","MARK1","MARK4","MBD1","MCL1","MDM2","MDM4","MEN1","MITF","MLL","MLL2","MLL3","MLLT10","MMP2","MN1","MRE11A","MTR","MTRR","MUC1","MUTYH","MYB","MYCL1","MYD88","MYH11","MYH9","NCOA1","NCOA2","NCOA4","NFKB1","NFKB2","NIN","NKX2-1","NLRP1","NOTCH4","NSD1","NUMA1","NUP214","NUP98","PAK3","PARP1","PAX3","PAX5","PAX7","PAX8","PBRM1","PBX1","PDE4DIP","PDGFB","PER1","PGAP3","PHOX2B","PIK3C2B","PIK3CD","PIK3CG","PIK3R2","PIM1","PKHD1","PLAG1","PLCG1","PLEKHG5","PML","PMS1","POT1","POU5F1","PPARG","PPP2R1A","PRDM1","PRKAR1A","PRKDC","PSIP1","PTGS2","PTPRD","PTPRT","RALGDS","RARA","RECQL4","REL","RHOH","RNASEL","RNF2","RNF213","RPS6KA2","RRM1","RUNX1T1","SAMD9","SBDS","SDHA","SDHB","SDHC","SOHD","SEPT9","SGK1","SH2D1A","SMAD2","SMAD4","SMUG1","SOCS1","SOX11","SOX2","SSX1","STK36","SUFU","SYK","SYNE1","TAF1","TAF1L","TAL1","TBX22","TCF12","TCF3","TCF7L1","TCF7L2","TCL1A","TET1","TFE3","TGFBR2","TGM7","THBS1","TIMP3","TLR4","TLX1","TNFAIP3","TNFRSF14","TNK2","TOP1","TPR","TRIM24","TRIM33","TRIP11","TRRAP","TSHR","UBR5","UGT1A1","USP9X","VHL","WAS","WHSC1","WRN","XPA","XPC","XPO1","XRCC2","ZNF384","ZNF521","ABL1","AKT1","AKT2","AKT3","ALK","AR","AXL","BRAF","CBL","CCND1","CDK4","CDK6","CSF1R","DDR2","EGFR","ERBB2","ERBB3","ERBB4","ERCC2","ESR1","EZH2","FGFR1","FGFR2","FGFR3","FGFR4","FLT3","GNA11","GNAQ","GNAS","HFN1A","HRAS","IDH1","IDH2","JAK2","KOR","KIT","KRAS","MAP2K1","MAP2K2","MAP2K4","MAPK1","MET","MPL","MTOR","MYC","MYCN","NFE2L2","NRAS","NTRK1","NTRK3","PDGFRA","PDGFRB","PIK3CA","PIK3CB","PTPN11","RAF1","RET","ROS1","SF3B1","SMO","SRC","ARID1A","ASXL1","ATM","ATR","ATRX","BAP1","CDK12","CDKN2A","CDKN2B","CEBPA","CHEK1","CHEK2","CREBBP","DNMT3A","FANCA","FANCD2","FBXW7","MLH1","MSH2","MSH6","NBN","NF1","NF2","NOTCH1","NOTCH2","NPM1","PALB2","PIK3R1","PMS2","PTCH1","PTEN","RADSO","RB1","RUNX1","SETD2","SMARCA4","SMARCB1","STK11","TET2","TP53","TSC1","TSC2","WT1"))
sum(OncoMineTMB)





