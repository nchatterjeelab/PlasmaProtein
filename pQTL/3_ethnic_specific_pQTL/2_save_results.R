
################################################################################
## check the examples where top SNPs is not in the other population
## (check their allele frequencies in 1000G)

################
# AA
rm(list=ls())

library(readr)
library(dplyr)

Black <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")

load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_AA.RData")
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
     file = "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_AA.RData")

## rare variants in AA
load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_AA.RData")
## two counts
(sum(mafin1000GEA$MAF <= 2/498/2)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3260734

Black <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")
m <- mafin1000GEA$MAF <= 2/498/2
ethnicspecific <- rep("FALSE", nrow(Black))
ethnicspecific[Notin1000G] <- "sentinel SNP is not in 1000Genome"
ethnicspecific[noSNPinWhite[!(noSNPinWhite %in% Notin1000G)][m]] <- "TRUE"
Black$ethnicspecific <- ethnicspecific
write_tsv(Black,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_2.0.txt")

marginEA$A1 == marginAA$A1
tmp <- marginEA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginAA$BETA) )
cor(tmp,marginAA$BETA)


################
# EA
rm(list=ls())

library(readr)
library(dplyr)

White <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")

load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/tmp_EA.RData")
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
     file = "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_EA.RData")

## rare variants in AA
load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/ethnic_specific_pQTL/results_EA.RData")
## two counts
(sum(mafin1000GAA$MAF <= 2/659/2)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09995002

White <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")
m <- mafin1000GAA$MAF <= 2/659/2
ethnicspecific <- rep("FALSE", nrow(White))
ethnicspecific[Notin1000G] <- "sentinel SNP is not in 1000Genome"
ethnicspecific[noSNPinBlack[!(noSNPinBlack %in% Notin1000G)][m]] <- "TRUE"
White$ethnicspecific <- ethnicspecific
write_tsv(White,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_2.0.txt")


marginEA$A1 == marginAA$A1
tmp <- marginAA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginEA$BETA) )
cor(tmp,marginEA$BETA)

