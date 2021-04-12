
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(dplyr)
library(stringr)

tissue <- dir("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS")
tissue <- tissue[str_detect(tissue,".P01.pos")]

n_gene <- integer()
for (i in 1:length(tissue)){
  tmp <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/", tissue[i]))
  n_gene <- c(n_gene,dim(tmp)[1])
}
names(n_gene) <- gsub("[.].*","",tissue)

saveRDS(n_gene, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/4_CompareTWASmodels/n_gene.rds")