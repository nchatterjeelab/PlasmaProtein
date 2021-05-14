
rm(list=ls())

library(readr)

for (chr in 1:22){
  if(chr==1){
    aask <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/50/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F)
  }else{
    aask <- rbind(aask, read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/50/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F))
  }
  print(chr)
}
for (chr in 1:22){
  if(chr==1){
    aric <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F)
  }else{
    aric <- rbind(aric, read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/permutation/chr",chr,"/permutations.txt"), stringsAsFactors = F))
  }
  print(chr)
}
aric_allsig <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/conditional/allsig.txt"),stringsAsFactors = F)
pGene_aric <- unique(aric_allsig$V1)
res_aask <- aask[aask$V1 %in% pgene_aric,]; dim(res)
# note: 1527 out of 1605 pGene_aric is in aask platform
summary(res_aask$V16)

res_aric <- aric[aric$V1 %in% pgene_aric,]; dim(res)
summary(res_aric$V16)



aask_allsig <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/50/conditional/allsig.txt"), stringsAsFactors = F)
pGene_aask <- unique(aask_allsig$V1) # 774
sum(pGene_aask %in% pGene_aric)
mean(pGene_aask %in% pGene_aric)




for (n_peer in (0:10)*10){

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
  #print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/allsig.txt"), col_names = F)
  print(paste0(n_peer,": ",length(unique(res$V1))))
}


#[1] "0: 508"
#[1] "10: 716"
#[1] "20: 740"
#[1] "30: 761"
#[1] "40: 770"
#[1] "50: 774" !!
#[1] "60: 751"
#[1] "70: 743"
#[1] "80: 721"
#[1] "90: 717"
#[1] "100: 712"


## clean results from conditonal analysis (exclude signals assigned with same SNPs)

cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
soma_tmp <- unique(cond$V1)
length(soma_tmp)
cond_cleaned <- tibble()
for (i in 1:length(soma_tmp)){
 cond_tmp <- cond[cond$V1 == soma_tmp[i],]
 if(nrow(cond_tmp)>1){
   rec_tmp <- integer()
   for (j in 2:nrow(cond_tmp)){
     for (r in 1:(j-1)){
       if(cond_tmp$V8[j] == cond_tmp$V8[r]){
         rec_tmp <- c(rec_tmp,j)
         break
       }
     }
   }
   if(length(rec_tmp) > 0){
     cond_tmp <- cond_tmp[-rec_tmp,]
   }
 }

 cond_cleaned <- rbind(cond_cleaned, cond_tmp)
 print(i)
}
write_tsv(cond_cleaned, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)




