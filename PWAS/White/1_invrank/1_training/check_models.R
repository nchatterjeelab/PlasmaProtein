
rm(list=ls())

library(readr)
library(stringr)

RDat <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_0.05/")
RDat <- RDat[str_detect(RDat,"RDat")]
seq <- unlist(lapply(str_split(RDat, "[.]wgt[.]RDat"), FUN=function (x){x[1]}))

for (i in 1:)