
rm(list=ls())

library(readr)
library(stringr)

n_peer <- 120

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
  print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), col_names = F)

length(unique(res$V1))

#0:  1515
#10: 1842
#20: 1894
#30: 1920
#40: 1883
#50: 1943
#60: 1943
#70: 1947
#80: 1952
#90: 1960
#100:1924
#110:1970
#120:1992 !!
#130:1948
#140:1975
#150:1974
#160:1939
#170:1942
#180:1956
#190:1924
#200:1935



rm(list=ls())

library(readr)
library(stringr)

n_peer <- 120

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt"))

  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/permutations_all.thresholds.txt"), col_names = F)


#rm(list=ls())
#
#library(readr)
#library(stringr)
#
#n_peer <- 120
#
#for (i in 1:22){
#  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt"))
#  tmp1 <- tmp[(tmp$V19==1) | (tmp$V15==1),]
#  if(i==1){
#    res <- tmp1
#  }else{
#    res <- rbind(res, tmp1)
#  }
#  print(i)
#}
#pGene.w <- unique(res$V1)
#
#res_new <- tibble()
#for (i in 1:length(pGene.w)){
#  tmp <- res[res$V1 == pGene.w[i],]
#  rank <- unique(tmp$V12)
#
#  for (j in 1:length(rank)){
#    if(sum(tmp$V12==rank[j])==1){
#      res_new <- rbind(res_new, tmp[tmp$V12==rank[j],])
#    }else{
#      res_new <- rbind(res_new, tmp[(tmp$V12==rank[j])&(tmp$V19==1),])
#    }
#  }
#  print(i)
#}
#


rm(list=ls())
library(readr)
n_peer <- 120
allsig <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"),stringsAsFactors = F)
allsig_sentinel <- allsig[allsig$V12 == 0,]
pGene.w <- unique(allsig$V1)

a <- integer()
for (i in 1:length(pGene.w)){
  tmp <- allsig[allsig$V1==pGene.w[i],]
  a[i] <- tmp$V12[1]
}

gene <- pGene.w[a==1]
for (i in 1:length(gene)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",allsig$V2[allsig$V1==gene[i]][1],"/conditional.txt"))
  allsig_sentinel <- rbind(allsig_sentinel, tmp[(tmp$V1==gene[i]) & (tmp$V12 == 0) & (tmp$V15 == 1) ,])
  print(i)
}

allsig_sentinel <- allsig_sentinel[order(allsig_sentinel$V3,decreasing=F),]
allsig_sentinel <- allsig_sentinel[order(allsig_sentinel$V2,decreasing=F),]
write_tsv(allsig_sentinel, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/topsig_corrected.txt"), col_names = F)

paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/conditional/topsig_corrected.txt")


#scp -r /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation  jzhan218@jhu.edu@login.marcc.jhu.edu:/home-2/jzhan218@jhu.edu/data/jzhan218/PWAS/Results/White/pQTL/
