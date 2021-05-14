
rm(list=ls())

# n_peer <- 60

for (i in 1:22){
  tmp <- dir(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i))
  if(length(tmp)==0)
  print(paste0("chr",i))
}


