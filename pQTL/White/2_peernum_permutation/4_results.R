
rm(list=ls())

library(readr)
library(stringr)

## number of pGenes for different number of  PEER factors

for (n_peer in (0:20)*10){

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

print(paste0(n_peer,": ",length(unique(res$V1))))

}

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
#120:1992 !! Best number of PEER factors
#130:1948
#140:1975
#150:1974
#160:1939
#170:1942
#180:1956
#190:1924
#200:1935


## significance threshold for all SOMAmers

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


## clean results from conditonal analysis (exclude signals assigned with same SNPs)

cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
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
write_tsv(cond_cleaned, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)

