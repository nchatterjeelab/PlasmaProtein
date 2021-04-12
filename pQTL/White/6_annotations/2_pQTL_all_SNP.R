
## white

rm(list=ls())

library(readr)
library(stringr)

seq <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/seqid_autosomal_withSNP.txt")

SNP <- character()
for (i in 1:length(seq)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/byseq/",seq[i],".bim"),stringsAsFactors = F)
  SNP <- unique(c(SNP, tmp$V2))
  print(i)
}
writeLines(SNP, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/all_cis_SNP.txt")

SNP <- unlist(lapply(SNP,FUN=function (x){substr(x,start=1, stop=nchar(x)-4)}))

a <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/GTex_results/WGS_Feature_overlap_collapsed_VEP_short_4torus_SNP.rds")
chr <- unlist(lapply(a, FUN=function (x){x[1]}))
pos <- unlist(lapply(a, FUN=function (x){x[2]}))
annotations_SNP <- paste0(chr,":", pos)

SNP_in <- SNP[SNP %in% annotations_SNP]
mean(SNP %in% annotations_SNP)



## black

rm(list=ls())

library(readr)
library(stringr)

seq <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/seqid_autosomal_withSNP.txt")

SNP <- character()
for (i in 1:length(seq)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/byseq/",seq[i],".bim"),stringsAsFactors = F)
  SNP <- unique(c(SNP, tmp$V2))
  print(i)
}
writeLines(SNP, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/all_cis_SNP.txt")
