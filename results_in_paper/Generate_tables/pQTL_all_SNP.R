
####################################
####################################
## pQTL summary

## white

rm(list=ls())

library(readr)
library(stringr)

seq <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/seqid_autosomal_withSNP.txt")
annota <- readr::read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
annota <- annota[,c(1,2,5,12,25:27)]
annota <- annota[annota$flag2==0,]
annota <- annota[!(is.na(annota$uniprot_id)),]

a <- annota[annota$seqid_in_sample %in% seq, ]

# nSNP <- integer()
# for (i in 1:22){
#   a <- system(paste0("wc -l /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr",i,".bim"), intern = TRUE)
#   a <- as.integer(gsub("[ ].*","", a))
#   nSNP[i] <- a
# }
#
# nSNP <- integer()
# for (i in 1:22){
#   a <- system(paste0("wc -l /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr",i,".bim"), intern = TRUE)
#   a <- as.integer(gsub("[ ].*","", a))
#   nSNP[i] <- a
# }
# sum(nSNP)

seq <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/seqid_autosomal_withSNP.txt")
SNP <- character()
for (i in 1:length(seq)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/byseq/",seq[i],".bim"))$V2
  SNP <- c(SNP, tmp)
  print(i)
}
SNP <- unique(SNP)
length(SNP) # [1] 6181856


seq <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/seqid_autosomal_withSNP.txt")
SNP <- character()
for (i in 1:length(seq)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/byseq/",seq[i],".bim"))$V2
  SNP <- c(SNP, tmp)
  print(i)
}
SNP <- unique(SNP)
length(SNP)

