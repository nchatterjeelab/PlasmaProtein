
################################################################################

rm(list=ls())

library(readr)
library(dplyr)

ethnic="Black"

PAV <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/3_PAV.txt"))
PAV <- PAV[,c(1,9,5,6,10)]

tab <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/1_pQTL_summary_cleaned_3.0_rsid.txt"))
tab <- left_join(tab, PAV, by=c("SOMAmer"="SeqID"))

write_tsv(tab, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/1_pQTL_summary_cleaned_4.0_rsid.txt"))

nrow(tab)

paste0(sum(tab$PAVLD=="R2 < 0.1"), " (",round(100*mean(tab$PAVLD=="R2 < 0.1"),1),"%)")

paste0(sum(tab$PAVLD=="0.1 <= R2 <= 0.9"), " (",round(100*mean(tab$PAVLD=="0.1 <= R2 <= 0.9"),1),"%)")

sum( (tab$PAVLD=="0.1 <= R2 <= 0.9") & (tab$sig) )

paste0(sum(tab$PAVLD=="R2 > 0.9"), " (",round(100*mean(tab$PAVLD=="R2 > 0.9"),1),"%)")

sum( (tab$PAVLD=="0.1 <= R2 <= 0.9") & (tab$sig) ) + sum(tab$PAVLD=="R2 > 0.9")
