
rm(list=ls())

library(readr)

for (chr in 1:22){
  if(chr==1){
    aask <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/40/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F)
  }else{
    aask <- rbind(aask, read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/40/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F))
  }
  print(chr)
}
for (chr in 1:22){
  if(chr==1){
    aric <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/80/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F)
  }else{
    aric <- rbind(aric, read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/80/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F))
  }
  print(chr)
}
aric_allsig <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/80/conditional/allsig.txt"),stringsAsFactors = F)
pGene_aric <- unique(aric_allsig$V1)
res_aask <- aask[aask$V1 %in% pGene_aric,]; dim(res_aask)
# note: 1541 out of 1618 pGene_aric is in aask platform


aask_allsig <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/40/conditional/allsig.txt"), stringsAsFactors = F)
pGene_aask <- unique(aask_allsig$V1) # 808
sum(pGene_aask %in% pGene_aric) # 761
mean(pGene_aask %in% pGene_aric) # 0.9418317


