
################################################################################
################################################################################

## pQTL 500kb

rm(list=ls())

library(readr)
library(stringr)


for (n_peer in (0:20)*10){

for (i in 1:22){
  tmp <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), col_names = F)

  print(paste0(n_peer,": ",length(unique(res$V1))))
}

#[1] "0: 1551"
#[1] "10: 1841"
#[1] "20: 1904"
#[1] "30: 1941"
#[1] "40: 1945"
#[1] "50: 1960"
#[1] "60: 1941"
#[1] "70: 1974"
#[1] "80: 1980"
#[1] "90: 2004"!
#[1] "100: 2004"
#[1] "110: 1996"
#[1] "120: 1989"
#[1] "130: 1966"
#[1] "140: 1970"
#[1] "150: 1975"
#[1] "160: 1953"
#[1] "170: 1971"
#[1] "180: 1974"
#[1] "190: 1989"
#[1] "200: 1982"



rm(list=ls())

library(readr)
library(stringr)


for (n_peer in (0:20)*10){

for (i in 1:22){
  tmp <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), col_names = F)

  print(paste0(n_peer,": ",length(unique(res$V1))))
}

#[1] "0: 1100"
#[1] "10: 1462"
#[1] "20: 1551"
#[1] "30: 1560"
#[1] "40: 1585"
#[1] "50: 1581"
#[1] "60: 1591"
#[1] "70: 1602"
#[1] "80: 1618"!
#[1] "90: 1593"
#[1] "100: 1589"
#[1] "110: 1601"
#[1] "120: 1588"
#[1] "130: 1544"
#[1] "140: 1549"
#[1] "150: 1551"
#[1] "160: 1561"
#[1] "170: 1554"
#[1] "180: 1564"
#[1] "190: 1533"
#[1] "200: 1519"


################################################################
################################################################

## significance threshold for all SOMAmers

rm(list=ls())

library(readr)
library(stringr)

n_peer <- 90

for (i in 1:22){
  tmp <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt"))

  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/permutations_all.thresholds.txt"), col_names = F)



rm(list=ls())

library(readr)
library(stringr)

n_peer <- 80

for (i in 1:22){
  tmp <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt"))

  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/permutations_all.thresholds.txt"), col_names = F)




########################################################################
########################################################################
## clean results from conditonal analysis (exclude signals assigned with same SNPs)

library(dplyr)

n_peer <- 90

cond <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
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
write_tsv(cond_cleaned, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)




library(dplyr)

n_peer <- 80

cond <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
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
write_tsv(cond_cleaned, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)





pGene.w <- unique(read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/90/conditional/allsig_cleaned.txt"))$V1)
pGene.b <- unique(read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/80/conditional/allsig_cleaned.txt"))$V1)
length(pGene.w) # 2004
length(pGene.b) # 1618
length(intersect(pGene.w,pGene.b)) # 1447
save(pGene.w, pGene.b, file="/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/genelist/genelist.rds")



pGene.w0 <- unique(read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/conditional/allsig_cleaned.txt"))$V1)
pGene.b0 <- unique(read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/conditional/allsig_cleaned.txt"))$V1)

length(setdiff(pGene.w,pGene.w0)) # 91
length(setdiff(pGene.b,pGene.b0)) # 90

newgene.w <- setdiff(pGene.w,pGene.w0)
newgene.b <- setdiff(pGene.b,pGene.b0)

save(newgene.w, newgene.b, file="/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/genelist/newgenelist.rds")


