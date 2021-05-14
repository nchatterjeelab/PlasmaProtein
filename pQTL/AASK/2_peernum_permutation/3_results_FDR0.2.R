
rm(list=ls())

library(readr)
#library(stringr)

for (n_peer in (0:10)*10){

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.2/",n_peer,"/conditional/chr",i,"/conditional.txt"))
  tmp1 <- tmp[tmp$V19==1,]
  if(i==1){
    res <- tmp1
  }else{
    res <- rbind(res, tmp1)
  }
  #print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.2/",n_peer,"/conditional/allsig.txt"), col_names = F)
  print(paste0(n_peer,": ",length(unique(res$V1))))
}

0.1
[1] "0: 580"
[1] "10: 806"
[1] "20: 834"
[1] "30: 841"
[1] "40: 856"
[1] "50: 872" !!
[1] "60: 869"
[1] "70: 859"
[1] "80: 798"
[1] "90: 798"
[1] "100: 812"

0.05
[1] "0: 508"
[1] "10: 716"
[1] "20: 740"
[1] "30: 761"
[1] "40: 770"
[1] "50: 774" !!
[1] "60: 751"
[1] "70: 743"
[1] "80: 721"
[1] "90: 717"
[1] "100: 712"



#0:  1088
#10: 1447
#20: 1507
#30: 1563
#40: 1593
#50: 1550
#60: 1574
#70: 1605 !! AFR
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
