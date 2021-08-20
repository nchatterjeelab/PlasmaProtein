
rm(list=ls())

library(readr)
library(plink2R)
library(bigreadr)
library(stringr)
library(dplyr)

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL")
White <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")
Black <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")

## EA
res_White <- tibble()
res_Black <- tibble()
noSNPinBlack <- integer() # index of top SNPs which are not in the other population (common MAF>1%)
for (i in 1:nrow(White)){
  seqid <- White$SOMAmer[i]
  topSNP <- White$TopSNP[i]
  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  res_White <- rbind(res_White,gwas[gwas$ID==topSNP,])

  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  if(topSNP %in% gwas$ID){
    res_Black <- rbind(res_Black,gwas[gwas$ID==topSNP,])
  }else{
    noSNPinBlack <- c(noSNPinBlack,i)
  }
  print(i)
}
save(noSNPinBlack,res_White,res_Black,
     file = "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_EA.RData")

## AA
res_White <- tibble()
res_Black <- tibble()
noSNPinWhite <- integer() # index of top SNPs which are not in the other population (common MAF>1%)
for (i in 1:nrow(Black)){
  seqid <- Black$SOMAmer[i]
  topSNP <- Black$TopSNP[i]
  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  res_Black <- rbind(res_Black,gwas[gwas$ID==topSNP,])

  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  if(topSNP %in% gwas$ID){
    res_White <- rbind(res_White,gwas[gwas$ID==topSNP,])
  }else{
    noSNPinWhite <- c(noSNPinWhite,i)
  }
  print(i)
}

save(noSNPinWhite,res_White,res_Black,
     file = "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_AA.RData")

