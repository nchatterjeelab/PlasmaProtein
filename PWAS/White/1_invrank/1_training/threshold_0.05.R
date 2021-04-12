
rm(list=ls())
library(stringr)

#system("cp -r /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_0.05")

all <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
hsq <- all[str_detect(all, "hsq")]; hsq <- gsub("[.].*","", hsq)
RDat <- all[str_detect(all, "RDat")]; RDat <-  gsub("[.].*","", RDat)
length(RDat)

dif <- setdiff(hsq, RDat)
new <- character()
for (i in 1:length(dif)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp/",dif[i], ".hsq"))
  if(tmp$V4 < 0.05){
    load(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_all_protein/",dif[i], ".wgt.RDat"))

    if(sum(wgt.matrix[,"enet"]!=0)>0){
      new <- c(new, dif[i])

      hsq <- as.numeric(tmp[1,2:3])
      hsq.pv <- tmp$V4
      save(cv.performance, hsq, hsq.pv, N.tot, snps, wgt.matrix,
           file = paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_0.05/",dif[i], ".wgt.RDat"))
    }

  }
}

length(new) #40


all <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_0.05/")
hsq <- all[str_detect(all, "hsq")]; hsq <- gsub("[.].*","", hsq); length(hsq)
RDat <- all[str_detect(all, "RDat")]; RDat <-  gsub("[.].*","", RDat); length(RDat)
