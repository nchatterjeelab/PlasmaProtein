
rm(list=ls())

library(readr)
library(stringr)

## number of pGenes for different number of  PEER factors

for (n_peer in (0:20)*10){

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
  print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), col_names = F)

print(paste0(n_peer,": ",length(unique(res$V1))))

}

#0:  1088
#10: 1447
#20: 1507
#30: 1563
#40: 1593
#50: 1550
#60: 1574
#70: 1605 !! Best number of PEER factor
#80: 1517
#90: 1566
#100:1340
#110:1564
#120:1320
#130:1480
#140:1474
#150:1534
#160:1523
#170: 1505
#180: 1492
#190: 1509
#200:1473


## significance threshold for all SOMAmers

rm(list=ls())

library(readr)
library(stringr)

n_peer <- 70 

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt"))

  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/permutations_all.thresholds.txt"), col_names = F)


## clean results from conditonal analysis (exclude signals assigned with same SNPs)

cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
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
write_tsv(cond_cleaned, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)




