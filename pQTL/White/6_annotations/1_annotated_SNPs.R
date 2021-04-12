

rm(list=ls())

library(readr)
library(dplyr)
library(stringr)


annotations <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/GTex_results/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt")
annotations <- annotations$SNP
a <- str_split(annotations, "_")
saveRDS(a, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/GTex_results/WGS_Feature_overlap_collapsed_VEP_short_4torus_SNP.rds")

chr <- unlist(lapply(a, FUN=function (x){x[1]}))
pos <- unlist(lapply(a, FUN=function (x){x[2]}))
annotations_SNP <- paste0(chr,":", pos)

allsig <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/120/conditional/allsig.txt", stringsAsFactors = F)
allsig <- paste0("chr",allsig$V9, ":", format(allsig$V10, trim = T, scientific = F))
mean(allsig %in% annotations_SNP)
allsig_in <- allsig[allsig %in% annotations_SNP]
