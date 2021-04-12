
rm(list=ls())


dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF")


for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N chr", i, "
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=30G
#$ -pe local 3
#$ -m e

module load R/3.6.1

/users/jzhang2/RESEARCH/tools/plink/plink2 \\
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr",i," \\
--freq \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",i,"

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/8_MAF/chr', i, '.sh'))


  print(i)
}


library(dplyr)
library(readr)
a <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/conditional/allsig.txt", stringsAsFactors = F)
res <- tibble()
for (chr in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",chr,".afreq"), stringsAsFactors = F)
  a.tmp <- a[a$V2 == chr,]
  res <- rbind(res,tmp[tmp$V2 %in% a.tmp$V8,])
  print(chr)
}
write_tsv(res, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/MAFsig.txt", col_names = F)



library(dplyr)
library(readr)
a <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/conditional/allsig.txt", stringsAsFactors = F)
res <- tibble()
for (chr in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/chr",chr,".afreq"), stringsAsFactors = F)
  a.tmp <- a[a$V2 == chr,]
  res <- rbind(res,tmp[tmp$V2 %in% a.tmp$V8,])
  print(chr)
}
write_tsv(res, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/MAFsig.txt", col_names = F)
