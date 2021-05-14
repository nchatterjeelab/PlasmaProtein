
rm(list=ls())

library(readr)
library(plink2R)
library(bigreadr)
library(stringr)
library(dplyr)

White <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")
Black <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")

## EA
res_White <- tibble()
res_Black <- tibble()
noSNPinBlack <- integer() # index of top SNPs which are not in the other population (common MAF>1%)
for (i in 1:nrow(White)){
  seqid <- White$SOMAmer[i]
  topSNP <- White$TopSNP[i]
  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  res_White <- rbind(res_White,gwas[gwas$ID==topSNP,])

  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  if(topSNP %in% gwas$ID){
    res_Black <- rbind(res_Black,gwas[gwas$ID==topSNP,])
  }else{
    noSNPinBlack <- c(noSNPinBlack,i)
  }
  print(i)
}
save(noSNPinBlack,res_White,res_Black,
     file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_EA.RData")

## AA
res_White <- tibble()
res_Black <- tibble()
noSNPinWhite <- integer() # index of top SNPs which are not in the other population (common MAF>1%)
for (i in 1:nrow(Black)){
  seqid <- Black$SOMAmer[i]
  topSNP <- Black$TopSNP[i]
  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  res_Black <- rbind(res_Black,gwas[gwas$ID==topSNP,])

  gwas <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",seqid,".PHENO1.glm.linear"),col_types = cols())
  if(topSNP %in% gwas$ID){
    res_White <- rbind(res_White,gwas[gwas$ID==topSNP,])
  }else{
    noSNPinWhite <- c(noSNPinWhite,i)
  }
  print(i)
}

save(noSNPinWhite,res_White,res_Black,
     file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_AA.RData")



################################################################################
## check the examples where top SNPs is not in the other population
## (check their allele frequencies in 1000G)

################
# AA
rm(list=ls())

library(readr)
library(dplyr)

Black <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")

load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_AA.RData")
annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')

SEQ <- Black$SOMAmer[noSNPinWhite]
CHR <- annota$chromosome_name[match(SEQ,annota$seqid_in_sample)]

Notin1000G <- integer()
reci <- integer()
res_EA <- tibble()
res_AA <- tibble()
for (chr in 1:22){
  m <- which(CHR == chr)
  afEA <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/EA_MAF_1000G/chr",chr,"_freq.afreq"), col_types = cols())
  afAA <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/AA_MAF_1000G/chr",chr,"_freq.afreq"), col_types = cols())

  for (j in 1:length(m)){
    i <- noSNPinWhite[m][j]
    reci <- c(reci,i)

    seqid <- Black$SOMAmer[i]
    topSNP <- Black$TopSNP[i]
    if( !(topSNP %in% afEA$ID) | !(topSNP %in% afAA$ID) ){
      Notin1000G <- c(Notin1000G,i)
      print("not in 1000G!")
      next
    }
    res_EA <- rbind(res_EA,afEA[afEA$ID == topSNP,])
    res_AA <- rbind(res_AA,afAA[afAA$ID == topSNP,])
    print(paste0("chr",chr,"--",i))
  }
}
tab1AA <- Black
marginEA <- res_White
marginAA <- res_Black[-noSNPinWhite,] ## top SNPs in the other population (common MAF>1%)
res_EA$MAF <- ifelse(res_EA$ALT_FREQS>0.5,1-res_EA$ALT_FREQS,res_EA$ALT_FREQS)
res_AA$MAF <- ifelse(res_AA$ALT_FREQS>0.5, 1-res_AA$ALT_FREQS, res_AA$ALT_FREQS)
## top SNPs' their allele frequencies/counts in 1000G
mafin1000GEA <- res_EA
mafin1000GAA <- res_AA
save(tab1AA, noSNPinWhite, marginEA, marginAA,
     reci, mafin1000GEA, mafin1000GAA, Notin1000G,
     file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_AA.RData")

## rare variants in EA
load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_AA.RData")
#(sum(mafin1000GEA$MAF < 0.01)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3641509
#(sum(mafin1000GEA$MAF < 0.005)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3427673
#(sum(mafin1000GEA$MAF < 0.002)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3
#(sum(mafin1000GEA$MAF == 0)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.236478
## two counts
(sum(mafin1000GEA$MAF <= 2/498/2)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3257862

Black <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")
m <- mafin1000GEA$MAF <= 2/659/2
ethnicspecific <- rep("FALSE", nrow(Black))
ethnicspecific[Notin1000G] <- "sentinel SNP is not in 1000Genome"
ethnicspecific[noSNPinWhite[!(noSNPinWhite %in% Notin1000G)][m]] <- "TRUE"
Black$ethnicspecific <- ethnicspecific
write_tsv(Black,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_2.0.txt")

marginEA$A1 == marginAA$A1
tmp <- marginEA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginAA$BETA) ) # 0.9495885
cor(tmp,marginAA$BETA) # 0.9521188


################
# EA
rm(list=ls())

library(readr)
library(dplyr)

White <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")

load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_EA.RData")
annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')

SEQ <- White$SOMAmer[noSNPinBlack]
CHR <- annota$chromosome_name[match(SEQ,annota$seqid_in_sample)]

Notin1000G <- integer()
reci <- integer()
res_EA <- tibble()
res_AA <- tibble()
for (chr in 1:22){
  m <- which(CHR == chr)
  afEA <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/EA_MAF_1000G/chr",chr,"_freq.afreq"), col_types = cols())
  afAA <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/AA_MAF_1000G/chr",chr,"_freq.afreq"), col_types = cols())

  for (j in 1:length(m)){
    i <- noSNPinBlack[m][j]
    reci <- c(reci,i)

    seqid <- White$SOMAmer[i]
    topSNP <- White$TopSNP[i]
    if( !(topSNP %in% afEA$ID) | !(topSNP %in% afAA$ID) ){
      Notin1000G <- c(Notin1000G,i)
      print("not in 1000G!")
      next
    }
    res_EA <- rbind(res_EA,afEA[afEA$ID == topSNP,])
    res_AA <- rbind(res_AA,afAA[afAA$ID == topSNP,])
    print(paste0("chr",chr,"--",i))
  }
}
tab1EA <- White
marginAA <- res_Black
marginEA <- res_White[-noSNPinBlack,] ## top SNPs in the other population (common MAF>1%)
res_AA$MAF <- ifelse(res_AA$ALT_FREQS>0.5, 1-res_AA$ALT_FREQS, res_AA$ALT_FREQS)
res_EA$MAF <- ifelse(res_EA$ALT_FREQS>0.5, 1-res_EA$ALT_FREQS, res_EA$ALT_FREQS)
## top SNPs' their allele frequencies/counts in 1000G
mafin1000GAA <- res_AA
mafin1000GEA <- res_EA
save(tab1EA, noSNPinBlack, marginAA, marginEA,
     reci, mafin1000GAA, mafin1000GEA, Notin1000G,
     file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_EA.RData")

## rare variants in AA
load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_EA.RData")
#(sum(mafin1000GAA$MAF < 0.01)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1242455
#(sum(mafin1000GAA$MAF < 0.005)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1222334
#(sum(mafin1000GAA$MAF < 0.002)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09507042
#(sum(mafin1000GAA$MAF == 0)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.03822938
## two counts
(sum(mafin1000GAA$MAF <= 2/659/2)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09507042

White <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")
m <- mafin1000GAA$MAF <= 2/659/2
ethnicspecific <- rep("FALSE", nrow(White))
ethnicspecific[Notin1000G] <- "sentinel SNP is not in 1000Genome"
ethnicspecific[noSNPinBlack[!(noSNPinBlack %in% Notin1000G)][m]] <- "TRUE"
White$ethnicspecific <- ethnicspecific
write_tsv(White,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_2.0.txt")


marginEA$A1 == marginAA$A1
tmp <- marginAA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginEA$BETA) ) # 0.9066282
cor(tmp,marginEA$BETA) # 0.9271907

############################################################

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
