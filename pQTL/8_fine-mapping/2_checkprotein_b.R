
################################################################
################################################################
################################################################

rm(list=ls())

#library(readr)

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt')
aveLD_b <- numeric()
for (i in 1:length(seqid)){
  gene <- seqid[i]
  bim <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/", gene, ".bim"))
  n.snp <- nrow(bim)
  LD <- readBin(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/LD/',gene, '.ld.bin'), what="numeric", size=4, n=(n.snp)^2)
  LD <- LD[LD!=1]
  aveLD_b[i] <- mean(LD^2)
  print(i)
}
saveRDS(aveLD_b, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/aveLD2_b.rds")



