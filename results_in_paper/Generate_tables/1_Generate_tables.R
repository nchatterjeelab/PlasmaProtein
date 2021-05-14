
####################################
####################################
## pQTL summary

## white

rm(list=ls())

library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 120

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.significant.txt"), stringsAsFactors = F)
  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), col_names = F)

res <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)

library(dplyr)
tab1 <- res[,c(1,2,3,6,8,10,17,16)]
colnames(tab1) <- c("SOMAmer", "Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval")
cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)
numcond <- as.data.frame(table(cond$V1), stringsAsFactors = F); colnames(numcond)[2] <- "NumCond"
tab1 <- inner_join(tab1, numcond, by=c("SOMAmer"="Var1"))
tab1 <- inner_join(tab1, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab1 <- tab1[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval","NumCond")]

write_tsv(tab1,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")


tab1 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned.txt")

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
write_tsv(tab1,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")


## Black

rm(list=ls())

library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 70

for (i in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.significant.txt"), stringsAsFactors = F)
  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
  print(i)
}
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), col_names = F)

res <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)
library(dplyr)
tab1 <- res[,c(1,2,3,6,8,10,17,16)]
colnames(tab1) <- c("SOMAmer", "Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval")
cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)
numcond <- as.data.frame(table(cond$V1), stringsAsFactors = F); colnames(numcond)[2] <- "NumCond"
tab1 <- inner_join(tab1, numcond, by=c("SOMAmer"="Var1"))
tab1 <- inner_join(tab1, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
                   by=c("SOMAmer"="seqid_in_sample"))
tab1 <- tab1[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","NumCisSNP","TopSNP","TopSNP_Pos38","Beta","Pval","NumCond")]

write_tsv(tab1,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")

# /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/SNPconvert/ARIC_GRCh38_ID.txt??


tab1 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")

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
write_tsv(tab1,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_1.0.txt")


####################################
####################################
## conditional report

## white

rm(list=ls())

library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 120

allpQTL <- read.table( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)
#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
#lost <- setdiff(unique(cond$V1), cond$V1[cond$V12==0]) # removed due to insignificance (this part is incorrect)
#length(lost)
#for (i in 1:length(lost)){
#  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",allpQTL$V2[allpQTL$V1==lost[i]],"/conditional.txt"),stringsAsFactors = F)
#  tmp <- tmp[(tmp$V1 == lost[i]) & (tmp$V12 == 1) ,]
#  tmp <- tmp[which.min(tmp$V17),]
#  a <- which(cond$V1 == lost[i])
#  cond <- rbind(cond[1:(a[1]-1),], tmp, cond[a[1]:(nrow(cond)),])
#  print(i)
#}
#write_tsv(cond, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_corrected.txt"), col_names = F)

#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
#soma_tmp <- unique(cond$V1)
#length(soma_tmp)
#cond_cleaned <- tibble()
#for (i in 1:length(soma_tmp)){
#  cond_tmp <- cond[cond$V1 == soma_tmp[i],]
#  if(nrow(cond_tmp)>1){
#    rec_tmp <- integer()
#    for (j in 2:nrow(cond_tmp)){
#      for (r in 1:(j-1)){
#        if(cond_tmp$V8[j] == cond_tmp$V8[r]){
#          rec_tmp <- c(rec_tmp,j)
#          break
#        }
#      }
#    }
#    if(length(rec_tmp) > 0){
#      cond_tmp <- cond_tmp[-rec_tmp,]
#    }
#  }
#
#  cond_cleaned <- rbind(cond_cleaned, cond_tmp)
#  print(i)
#}
#write_tsv(cond_cleaned, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)


#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)
#
#library(dplyr)
#tab2 <- cond[,c(1,2,3,8,10,12,18,17)]
#colnames(tab2) <- c("SOMAmer", "Chr","TSS","SNP","SNP_Pos38","Rank","Beta","Pval")
#tmp <- allpQTL[,c(1,8,10)]; colnames(tmp) <- c("SOMAmer","TopSNP","TopSNP_Pos38")
#tab2 <- inner_join( tab2, tmp, by="SOMAmer")
#tab2 <- tab2[,c(1:3,9:10,4:8)]
#tab2 <- inner_join(tab2, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
#                   by=c("SOMAmer"="seqid_in_sample"))
#tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","Beta","Pval")]

#SOMA <- unique(tab2$SOMAmer)
#library( snpStats)
#length(SOMA)
#R2 <- numeric()
#for (i in 1:length(SOMA)){
#  tmp <- tab2[tab2$SOMAmer==SOMA[i],]
#  for (j in 1:nrow(tmp)){
#    if(tmp$TopSNP[1]  == tmp$SNP[j]){
#      R2 <- c(R2, 1)
#    }else{
#      a <- read.plink(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr",tmp$Chr[1]),select.snps=c(tmp$TopSNP[1], tmp$SNP[j]))
#      a <- as(a$genotypes, Class="numeric")
#      if(sum(is.na(a))!=0){
#        a[,1][is.na(a[,1])] <- mean(a[,1], na.rm = T)
#        a[,2][is.na(a[,2])] <- mean(a[,2], na.rm = T)
#      }
#      R2 <- c(R2, cor(a[,1],a[,2])^2)
#    }
#  }
#  print(i)
#}
#tab2 <- cbind(tab2, R2)
#!!! see 2_tab2 folder
tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/tab2/tab2_R2.txt")

#nominal <- list()
#for (chr in 1:22){
#  nominal[[chr]] <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/chr",chr,".txt"), stringsAsFactors = F)
#  print(chr)
#}
#saveRDS(nominal, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/all.rds")
nominal <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/nominal/all.rds")

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
write_tsv(tab2,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis.txt")

tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis.txt")
tab2 <- tab2[tab2$Rank !=0,]
length(unique(tab2$SOMAmer))

tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis.txt")
mean(abs(tab2$TopSNP_Pos38-tab2$SNP_Pos38)<200000)

tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis.txt")
sum(tab2$Rank==5)

tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis.txt")
sum(abs(tab2$Beta>0.5))


## black

rm(list=ls())

library(readr)
library(stringr)
library(dplyr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 70

allpQTL <- read.table( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)
#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
#lost <- setdiff(unique(cond$V1), cond$V1[cond$V12==0])
#length(lost)
#for (i in 1:length(lost)){
#  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/chr",allpQTL$V2[allpQTL$V1==lost[i]],"/conditional.txt"),stringsAsFactors = F)
#  tmp <- tmp[(tmp$V1 == lost[i]) & (tmp$V12 == 1) ,]
#  tmp <- tmp[which.min(tmp$V17),]
#  a <- which(cond$V1 == lost[i])
#  cond <- rbind(cond[1:(a[1]-1),], tmp, cond[a[1]:(nrow(cond)),])
#  print(i)
#}
#write_tsv(cond, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_corrected.txt"), col_names = F)
#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"), stringsAsFactors = F)
#soma_tmp <- unique(cond$V1)
#length(soma_tmp)
#cond_cleaned <- tibble()
#for (i in 1:length(soma_tmp)){
#  cond_tmp <- cond[cond$V1 == soma_tmp[i],]
#  if(nrow(cond_tmp)>1){
#    rec_tmp <- integer()
#    for (j in 2:nrow(cond_tmp)){
#      for (r in 1:(j-1)){
#        if(cond_tmp$V8[j] == cond_tmp$V8[r]){
#          rec_tmp <- c(rec_tmp,j)
#          break
#        }
#      }
#    }
#    if(length(rec_tmp) > 0){
#      cond_tmp <- cond_tmp[-rec_tmp,]
#    }
#  }
#
#  cond_cleaned <- rbind(cond_cleaned, cond_tmp)
#  print(i)
#}
#write_tsv(cond_cleaned, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), col_names = F)


#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_corrected.txt"), stringsAsFactors = F)

#library(dplyr)
#tab2 <- cond[,c(1,2,3,8,10,12,18,17)]
#colnames(tab2) <- c("SOMAmer", "Chr","TSS","SNP","SNP_Pos38","Rank","Beta","Pval")
#tmp <- allpQTL[,c(1,8,10)]; colnames(tmp) <- c("SOMAmer","TopSNP","TopSNP_Pos38")
#tab2 <- inner_join( tab2, tmp, by="SOMAmer")
#tab2 <- tab2[,c(1:3,9:10,4:8)]
#tab2 <- inner_join(tab2, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
#                   by=c("SOMAmer"="seqid_in_sample"))
#tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","Beta","Pval")]
#
#SOMA <- unique(tab2$SOMAmer)
#library( snpStats)
#length(SOMA)
#R2 <- numeric()
#for (i in 1:length(SOMA)){
#  tmp <- tab2[tab2$SOMAmer==SOMA[i],]
#  for (j in 1:nrow(tmp)){
#    if(tmp$TopSNP[1]  == tmp$SNP[j]){
#      R2 <- c(R2, 1)
#    }else{
#      a <- read.plink(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr",tmp$Chr[1]),select.snps=c(tmp$TopSNP[1], tmp$SNP[j]))
#      a <- as(a$genotypes, Class="numeric")
#      if(sum(is.na(a))!=0){
#        a[,1][is.na(a[,1])] <- mean(a[,1], na.rm = T)
#        a[,2][is.na(a[,2])] <- mean(a[,2], na.rm = T)
#      }
#      R2 <- c(R2, cor(a[,1],a[,2])^2)
#    }
#  }
#  print(i)
#}
#tab2 <- cbind(tab2, R2)
# see 2_tab2 folder
tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/tab2/tab2_R2.txt")

#nominal <- list()
#for (chr in 1:22){
#  nominal[[chr]] <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/chr",chr,".txt"), stringsAsFactors = F)
#  print(chr)
#}
#saveRDS(nominal, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/all.rds")
nominal <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/nominal/all.rds")

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
write_tsv(tab2,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis.txt")

library(readr)
tab2 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis.txt")
a <- tab2$SNP
b <- bigreadr::fread2("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/rsid/GRCh38_dbSNP151_rsid_final_USE.txt")

