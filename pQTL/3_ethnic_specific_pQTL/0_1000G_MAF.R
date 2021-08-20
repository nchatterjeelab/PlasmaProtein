
rm(list=ls())

library(readr)
library(plink2R)
library(bigreadr)
library(stringr)
library(dplyr)

lookup <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")

for (chr in 1:22){
  af <- fread2(paste0("/dcl01/chatterj/data/jzhang2/1000G/GRCh38/EUR/chr",chr,"_freq.afreq"))
  af <- af[af$ID %in% lookup$rsid,]
  af$ID <- lookup$SNPid[match(af$ID,lookup$rsid)]
  write_tsv(af,paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/EA_MAF_1000G/chr",chr,"_freq.afreq"))
  print(chr)
}

for (chr in 1:22){
  af <- fread2(paste0("/dcl01/chatterj/data/jzhang2/1000G/GRCh38/AFR/chr",chr,"_freq.afreq"))
  af <- af[af$ID %in% lookup$rsid,]
  af$ID <- lookup$SNPid[match(af$ID,lookup$rsid)]
  write_tsv(af,paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/AA_MAF_1000G/chr",chr,"_freq.afreq"))
  print(chr)
}
