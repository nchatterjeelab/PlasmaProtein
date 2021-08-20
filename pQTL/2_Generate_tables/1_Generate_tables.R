
####################################
####################################
## pQTL summary

## white

rm(list=ls())

library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 90

for (i in 1:22){
  tmp <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.significant.txt"), stringsAsFactors = F)
  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), col_names = F)

res <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)

library(dplyr)
tab1 <- res[,c(1,2,3,6,8,10,17,16)]
colnames(tab1) <- c("SOMAmer", "Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval")
cond <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)
numcond <- as.data.frame(table(cond$V1), stringsAsFactors = F); colnames(numcond)[2] <- "NumCond"
tab1 <- inner_join(tab1, numcond, by=c("SOMAmer"="Var1"))
tab1 <- inner_join(tab1, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab1 <- tab1[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval","NumCond")]

write_tsv(tab1,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")


tab1 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")

ALT <- character()
AF <- numeric()
i=0
for (chr in 1:22){
  info <- read_tsv(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/info/chr",chr,".info.txt"), col_types = cols())
  maf <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/chr",chr,".afreq"), col_types = cols())
  m <- which(tab1$Chr == chr)
  for (j in 1:length(m)){
    i <- i+1; print(i)

    ALT[i] <- info$ALT[info$ID == tab1$TopSNP[m[j]]]
    tmp <- maf[maf$ID == tab1$TopSNP[m[j]],]
    if(tmp$ALT == ALT[i]){
      AF[i] <- tmp$ALT_FREQS
    }else{
      AF[i] <- 1-tmp$ALT_FREQS
    }
  }
}

tab1$A1 <- ALT
tab1$A1_AF <- AF
write_tsv(tab1,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")



## Black

rm(list=ls())

library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 80

for (i in 1:22){
  tmp <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.significant.txt"), stringsAsFactors = F)
  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), col_names = F)

res <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)

library(dplyr)
tab1 <- res[,c(1,2,3,6,8,10,17,16)]
colnames(tab1) <- c("SOMAmer", "Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval")
cond <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)
numcond <- as.data.frame(table(cond$V1), stringsAsFactors = F); colnames(numcond)[2] <- "NumCond"
tab1 <- inner_join(tab1, numcond, by=c("SOMAmer"="Var1"))
tab1 <- inner_join(tab1, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab1 <- tab1[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval","NumCond")]

write_tsv(tab1,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")


tab1 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")

ALT <- character()
AF <- numeric()
i=0
for (chr in 1:22){
  info <- read_tsv(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/info/chr",chr,".info.txt"), col_types = cols())
  maf <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",chr,".afreq"), col_types = cols())
  m <- which(tab1$Chr == chr)
  for (j in 1:length(m)){
    i <- i+1; print(i)

    ALT[i] <- info$ALT[info$ID == tab1$TopSNP[m[j]]]
    tmp <- maf[maf$ID == tab1$TopSNP[m[j]],]
    if(tmp$ALT == ALT[i]){
      AF[i] <- tmp$ALT_FREQS
    }else{
      AF[i] <- 1-tmp$ALT_FREQS
    }
  }
}

tab1$A1 <- ALT
tab1$A1_AF <- AF
write_tsv(tab1,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")


### 1_pQTL_summary_cleaned_2.0.txt please see  3_ethnic_specific_pQTL!!!

####################################
####################################
## conditional report

## white

rm(list=ls())

library(readr)
library(stringr)
n_peer <- 90
annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
allpQTL <- read.table( paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)
cond <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)

library(dplyr)
tab2 <- cond[,c(1,2,3,8,10,12,18,17)]
colnames(tab2) <- c("SOMAmer", "Chr","TSS","SNP","SNP_Pos38","Rank","Beta","Pval")
tmp <- allpQTL[,c(1,8,10)]; colnames(tmp) <- c("SOMAmer","TopSNP","TopSNP_Pos38")
tab2 <- inner_join( tab2, tmp, by="SOMAmer")
tab2 <- tab2[,c(1:3,9:10,4:8)]
tab2 <- inner_join(tab2, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","Beta","Pval")]
write_tsv(tab2,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/tab2.txt")
# R2 colum: see 2_tab2 folder!!!
tab2 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/tab2/tab2_R2.txt")

#nominal <- list()
#for (chr in 1:22){
#  nominal[[chr]] <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/chr",chr,".txt"), stringsAsFactors = F)
#  print(chr)
#}
#saveRDS(nominal, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/all.rds")
nominal <- readRDS("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/all.rds")

Marg_Beta <- numeric()
Marg_p <- numeric()
for (i in 1:nrow(tab2)){
  tmp <- nominal[[tab2$Chr[i]]]
  tmp <- tmp[tmp$V1 == tab2$SOMAmer[i],]
  tmp <- tmp[tmp$V8 == tab2$SNP[i],]
  Marg_Beta[i] <- tmp$V13
  Marg_p[i] <- tmp$V12
  print(i)
}
tab2 <- cbind(tab2, Marg_Beta, Marg_p)
tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","R2","Beta","Pval","Marg_Beta", "Marg_p")]
write_tsv(tab2,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis.txt")


library(readr)
tab2 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis.txt")

ALT <- character()
AF <- numeric()
i=0
for (chr in 1:22){
  info <- read_tsv(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/info/chr",chr,".info.txt"), col_types = cols())
  maf <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/chr",chr,".afreq"), col_types = cols())
  m <- which(tab2$Chr == chr)
  for (j in 1:length(m)){
    i <- i+1; print(i)

    ALT[i] <- info$ALT[info$ID == tab2$SNP[m[j]]]
    tmp <- maf[maf$ID == tab2$SNP[m[j]],]
    if(tmp$ALT == ALT[i]){
      AF[i] <- tmp$ALT_FREQS
    }else{
      AF[i] <- 1-tmp$ALT_FREQS
    }
  }
}
tab2$A1 <- ALT
tab2$A1_AF <- AF
write_tsv(tab2,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis_1.0.txt")


## black


rm(list=ls())

library(readr)
library(stringr)
n_peer <- 80
annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
allpQTL <- read.table( paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)
cond <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)

library(dplyr)
tab2 <- cond[,c(1,2,3,8,10,12,18,17)]
colnames(tab2) <- c("SOMAmer", "Chr","TSS","SNP","SNP_Pos38","Rank","Beta","Pval")
tmp <- allpQTL[,c(1,8,10)]; colnames(tmp) <- c("SOMAmer","TopSNP","TopSNP_Pos38")
tab2 <- inner_join( tab2, tmp, by="SOMAmer")
tab2 <- tab2[,c(1:3,9:10,4:8)]
tab2 <- inner_join(tab2, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","Beta","Pval")]
write_tsv(tab2,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/tab2.txt")
# R2 colum: see 2_tab2 folder!!!
tab2 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/tab2/tab2_R2.txt")


#nominal <- list()
#for (chr in 1:22){
#  nominal[[chr]] <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/chr",chr,".txt"), stringsAsFactors = F)
#  print(chr)
#}
#saveRDS(nominal, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/all.rds")
nominal <- readRDS("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/all.rds")

Marg_Beta <- numeric()
Marg_p <- numeric()
for (i in 1:nrow(tab2)){
  tmp <- nominal[[tab2$Chr[i]]]
  tmp <- tmp[tmp$V1 == tab2$SOMAmer[i],]
  tmp <- tmp[tmp$V8 == tab2$SNP[i],]
  Marg_Beta[i] <- tmp$V13
  Marg_p[i] <- tmp$V12
  print(i)
}
tab2 <- cbind(tab2, Marg_Beta, Marg_p)
tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","R2","Beta","Pval","Marg_Beta", "Marg_p")]
write_tsv(tab2,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis.txt")


library(readr)
tab2 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis.txt")

ALT <- character()
AF <- numeric()
i=0
for (chr in 1:22){
  info <- read_tsv(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/info/chr",chr,".info.txt"), col_types = cols())
  maf <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",chr,".afreq"), col_types = cols())
  m <- which(tab2$Chr == chr)
  for (j in 1:length(m)){
    i <- i+1; print(i)

    ALT[i] <- info$ALT[info$ID == tab2$SNP[m[j]]]
    tmp <- maf[maf$ID == tab2$SNP[m[j]],]
    if(tmp$ALT == ALT[i]){
      AF[i] <- tmp$ALT_FREQS
    }else{
      AF[i] <- 1-tmp$ALT_FREQS
    }
  }
}

tab2$A1 <- ALT
tab2$A1_AF <- AF
write_tsv(tab2,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis_1.0.txt")

