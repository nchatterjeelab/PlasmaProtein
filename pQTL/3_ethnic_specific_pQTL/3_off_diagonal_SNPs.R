
############################################################
############################################################
### check off-diagonal SNPs ###

#######################
## EA
load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_EA.RData")

tmp <- marginAA$BETA; tmp[marginAA$A1 != marginEA$A1] <- - tmp[marginAA$A1 != marginEA$A1]
mean( sign(tmp) == sign(marginEA$BETA) ) # 0.9495885
cor(tmp,marginEA$BETA) # 0.9521188
df.EA <- data.frame(Beta_AA=tmp,Beta_EA=marginEA$BETA,ID=marginEA$ID,stringsAsFactors = F)

library(readr)
library(dplyr)
allID <- marginAA$ID
chr <- gsub("chr","",unlist(lapply( stringr::str_split(allID,":"), FUN = function (x){x[1]})))
uni_chr <- unique(chr)
MAF_AA <- tibble()
for (i in 1:length(uni_chr)){
  tmp <- suppressMessages(read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/chr",uni_chr[i],".afreq")))
  m <- which(chr==uni_chr[i])
  MAF_AA <- rbind(MAF_AA, tmp[match(allID[m],tmp$ID),])
  print(i)
}
MAF_EA <- tibble()
for (i in 1:length(uni_chr)){
  tmp <- suppressMessages(read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",uni_chr[i],".afreq")))
  m <- which(chr==uni_chr[i])
  MAF_EA <- rbind(MAF_EA, tmp[match(allID[m],tmp$ID),])
  print(i)
}
MAF_AA$ALT==MAF_EA$ALT
MAF_AA$ALT_FREQS/MAF_EA$ALT_FREQS
#> summary(MAF_AA$ALT_FREQS/MAF_EA$ALT_FREQS)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03362 0.34991 0.91709 1.54944 1.80245 6.47576
save(df.EA,allID,MAF_AA,MAF_EA,file="/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/AF_check_EA.RData")

df.EA[(df.EA$Beta_EA > 0.5)&(df.EA$Beta_AA < -0.5),]
MAF_EA[MAF_EA$ID=="chr16:72045758:C:T",]
MAF_AA[MAF_AA$ID=="chr16:72045758:C:T",]

df.EA[(df.EA$Beta_EA < -1)&(df.EA$Beta_AA > 0),]
MAF_EA[MAF_EA$ID=="chr17:35980421:T:C",]
MAF_AA[MAF_AA$ID=="chr17:35980421:T:C",]

df.EA[(df.EA$Beta_EA < -0.5)&(df.EA$Beta_AA > 0.3),]
MAF_EA[MAF_EA$ID=="chr9:127878806:G:T",]
MAF_AA[MAF_AA$ID=="chr9:127878806:G:T",]

df.EA[(df.EA$Beta_EA > 1) & (df.EA$Beta_AA < 0.3),]
MAF_EA[MAF_EA$ID=="chr2:162209363:G:A",]
MAF_AA[MAF_AA$ID=="chr2:162209363:G:A",]


#######################
## AA
load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_AA.RData")

tmp <- marginEA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginAA$BETA) ) # 0.9495885
cor(tmp,marginAA$BETA) # 0.9521188
df.AA <- data.frame(Beta_EA=tmp,Beta_AA=marginAA$BETA,ID=marginAA$ID,stringsAsFactors = F)

library(readr)
library(dplyr)
allID <- marginEA$ID
chr <- gsub("chr","",unlist(lapply( stringr::str_split(allID,":"), FUN = function (x){x[1]})))
uni_chr <- unique(chr)
MAF_EA <- tibble()
for (i in 1:length(uni_chr)){
  tmp <- suppressMessages(read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/chr",uni_chr[i],".afreq")))
  m <- which(chr==uni_chr[i])
  MAF_EA <- rbind(MAF_EA, tmp[match(allID[m],tmp$ID),])
  print(i)
}
MAF_AA <- tibble()
for (i in 1:length(uni_chr)){
  tmp <- suppressMessages(read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",uni_chr[i],".afreq")))
  m <- which(chr==uni_chr[i])
  MAF_AA <- rbind(MAF_AA, tmp[match(allID[m],tmp$ID),])
  print(i)
}
MAF_EA$ALT==MAF_AA$ALT
MAF_EA$ALT_FREQS/MAF_AA$ALT_FREQS
#> summary(MAF_EA$ALT_FREQS/MAF_AA$ALT_FREQS)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03362 0.34991 0.91709 1.54944 1.80245 6.47576
save(df.AA,allID,MAF_EA,MAF_AA,
     file="/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/AF_check_AA.RData")

df.AA[(df.AA$Beta_EA<0.8)&(df.AA$Beta_AA>1),]
MAF_EA[MAF_EA$ID=="chr14:94509525:C:T",]
MAF_AA[MAF_AA$ID=="chr14:94509525:C:T",]

df.AA[(df.AA$Beta_EA < -0.2) & (df.AA$Beta_AA > 0.4),]
MAF_EA[MAF_EA$ID=="chr1:153346340:G:T",]
MAF_AA[MAF_AA$ID=="chr1:153346340:G:T",]

df.AA[(df.AA$Beta_EA > 0) & (df.AA$Beta_AA < -0.5),]
MAF_EA[MAF_EA$ID=="chr19:51627866:A:G",]
MAF_AA[MAF_AA$ID=="chr19:51627866:A:G",]
MAF_EA[MAF_EA$ID=="chr6:31852950:C:T",]
MAF_AA[MAF_AA$ID=="chr6:31852950:C:T",]

df.AA[(df.AA$Beta_EA > 1) & (df.AA$Beta_AA < 0.5),]
MAF_EA[MAF_EA$ID=="chr2:162224922:C:T",]
MAF_AA[MAF_AA$ID=="chr2:162224922:C:T",]
