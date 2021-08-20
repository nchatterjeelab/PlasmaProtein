
################################################################################

rm(list=ls())

library(readr)
library(stringr)

ethnic="Black"


for (i in 1:22){
  tmp <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/allsig.txt"), col_names = F)


res <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/allsig.txt"), stringsAsFactors = F)
length(unique(res$V1)) # 1801
a <- unique(res$V1)

original0.05 <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/1_pQTL_summary_cleaned_2.0_rsid.txt"))
backup0.01 <- cbind(original0.05, backup0.01=original0.05$SOMAmer %in% a)
write_tsv(backup0.01, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/1_pQTL_summary_cleaned_3.0_rsid.txt"))


