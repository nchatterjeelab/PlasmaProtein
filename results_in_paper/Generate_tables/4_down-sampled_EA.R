
rm(list=ls())

library(readr)
library(dplyr)
library(stringr)


###############################################################
## match sample size

cond.b <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/conditional/allsig.txt",stringsAsFactors = F)
cond.w <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/120/conditional/allsig.txt", stringsAsFactors = F)

incommon <- logical()
p <- numeric()
for (chr in 1:22){
  nomi.w <- bigreadr::fread2(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/9_White_matchNblack/all_sample_peers/nominal/chr",chr,".txt"))
  cond.b.tmp <- cond.b[cond.b$V2==chr,]
  incommon <- c(incommon, cond.b.tmp$V8 %in% nomi.w$V8)

  cond.b.tmp <- cond.b.tmp[cond.b.tmp$V8 %in% nomi.w$V8, ]
  for (i in 1:nrow(cond.b.tmp)){
    gene <- cond.b.tmp$V1[i]
    SNP <- cond.b.tmp$V8[i]
    nomi.tmp <- nomi.w[nomi.w$V1==gene,]
    nomi.tmp <- nomi.tmp[nomi.tmp$V8==SNP,]
    p <- c(p,nomi.tmp$V12)
  }
  print(chr)
}
save(incommon, p, file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/9_White_matchNblack/all_sample_peers/nominal/binw.rdata")

length(incommon) #4245
sum(incommon) # 2348
mean(incommon) # 0.5531213

length(p) # 2348
sum(p<0.001) # 1741
mean(p<0.001) # 0.7414821


incommon <- logical()
p <- numeric()
for (chr in 1:22){
  nomi.b <- bigreadr::fread2(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/chr",chr,".txt"))
  cond.w.tmp <- cond.w[cond.w$V2==chr,]
  incommon <- c(incommon, cond.w.tmp$V8 %in% nomi.b$V8)

  cond.w.tmp <- cond.w.tmp[cond.w.tmp$V8 %in% nomi.b$V8, ]
  for (i in 1:nrow(cond.w.tmp)){
    gene <- cond.w.tmp$V1[i]
    SNP <- cond.w.tmp$V8[i]
    nomi.tmp <- nomi.b[nomi.b$V1==gene,]
    nomi.tmp <- nomi.tmp[nomi.tmp$V8==SNP,]
    p <- c(p,nomi.tmp$V12)
  }
  print(chr)
}
save(incommon, p, file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/9_White_matchNblack/all_sample_peers/nominal/winb.rdata")

length(incommon) # 3005
sum(incommon) # 2492
mean(incommon) # 0.8292845

length(p) # 2492
sum(p<0.001) # 1488
mean(p<0.001) # 0.5971108


a <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/120/conditional/allsig.txt", stringsAsFactors = F)
b <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/conditional/allsig.txt", stringsAsFactors = F)
dim(a) # 3005
dim(b) # 7638

a1 <- paste0(a$V1,"-",a$V8)
b1 <- paste0(b$V1,"-",b$V8)


###

## replicate

rm(list=ls())

library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 120

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/120/permutation/chr",i,"/permutations_all.significant.txt"), stringsAsFactors = F)
  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}

library(dplyr)
tab1 <- res[,c(1,2,3,6,8,10,17,16)]
colnames(tab1) <- c("SOMAmer", "Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval")
cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/120/conditional/allsig.txt"), stringsAsFactors = F)
numcond <- as.data.frame(table(cond$V1), stringsAsFactors = F); colnames(numcond)[2] <- "NumCond"
tab1 <- inner_join(tab1, numcond, by=c("SOMAmer"="Var1"))
tab1 <- inner_join(tab1, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab1 <- tab1[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval","NumCond")]

p <- list()
for (chr in 1:22){
  tab1.tmp <- tab1[tab1$Chr==chr,]
  nomi.w <- bigreadr::fread2(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/chr",chr,".txt"))
  p[[chr]] <- dplyr::inner_join(tab1.tmp,nomi.w, by=c("SOMAmer"="V1", "TopSNP"="V8"))$V12
  print(chr)
}
p <- unlist(p)
tab1$P_allsample <- p
readr::write_tsv(tab1, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/9_White_matchNblack/all_sample_peers/tab1.txt")
mean(p<0.05/nrow(tab1)) # 0.925671
p_adj <- p.adjust(p, method = "fdr")
mean(p_adj<0.05) # 0.9662767




tab1 <- readr::read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary.txt")
p <- list()
for (chr in 1:22){
  tab1.tmp <- tab1[tab1$Chr==chr,]
  nomi.w <- bigreadr::fread2(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/9_White_matchNblack/all_sample_peers/nominal/chr",chr,".txt"))
  p[[chr]] <- dplyr::inner_join(tab1.tmp,nomi.w, by=c("SOMAmer"="V1", "TopSNP"="V8"))$V12
  print(chr)
}
p <- unlist(p)
mean(p<0.05/nrow(tab1)) # 0.6511044
p_adj <- p.adjust(p, method = "fdr")
mean(p_adj<0.05) # 0.8970884
