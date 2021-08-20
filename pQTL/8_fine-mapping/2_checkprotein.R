
################################################################
################################################################
################################################################

rm(list=ls())

#library(readr)

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/seqid_autosomal_withSNP.txt')
aveLD_w <- numeric()
for (i in 1:length(seqid)){
  gene <- seqid[i]
  bim <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/", gene, ".bim"))
  n.snp <- nrow(bim)
  LD <- readBin(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/LD/',gene, '.ld.bin'), what="numeric", size=4, n=(n.snp)^2)
  LD <- LD[LD!=1]
  aveLD_w[i] <- mean(LD^2)
  print(i)
}
saveRDS(aveLD_w, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/aveLD2_w.rds")



aveLD_w <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/aveLD2_w.rds")
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/seqid_autosomal_withSNP.txt')
names(aveLD_w) <- seqid

aveLD_b <-  readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/aveLD2_b.rds")
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt')
names(aveLD_b) <- seqid

pQTL.w <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/conditional/allsig.txt')
pGene.w <- unique(pQTL.w$V1)

pQTL.b <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/conditional/allsig.txt')
pGene.b <- unique(pQTL.b$V1)

seqid <- intersect(pGene.w,pGene.b)

df <- data.frame(seqid=seqid, aveLD_w=aveLD_w[seqid], aveLD_b=aveLD_b[seqid])
df$rank_w <- rank(df$aveLD_w)
df$rank_b <- rank(df$aveLD_b)
df$rank_diff <- df$rank_w - df$rank_b
rownames(df) <- NULL
df <- df[order(df$rank_w,decreasing=T),]; head(df)
df <- df[order(df$rank_diff,decreasing=T),]; head(df)

SeqId_14205_6,
SeqId_15596_7,
SeqId_14203_3,
SeqId_18244_1,
SeqId_9380_2

SeqId_8007_19,
SeqId_18343_10,
SeqId_16809_1,
SeqId_6919_3,
SeqId_13720_95



