
rm(list=ls())

#library(readr)
#library(stringr)

n_peer <- 60

for (n_peer in (1:10)*10){
for (i in 1:22){
  tmp <- dir(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/chr",i))
  if(length(tmp)==0)
  print(paste0(n_peer,"chr",i))
}
}
#write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), col_names = F)



