
####################################
####################################
## conditional report

## white

rm(list=ls())

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

i  <- as.integer(i)


library(readr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')

n_peer <- 120

allpQTL <- read.table( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/allpQTL.txt"), stringsAsFactors = F)
#cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_corrected.txt"), stringsAsFactors = F)
cond <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/allsig_cleaned.txt"), stringsAsFactors = F)

library(dplyr)
#tab2 <- cond[,c(1,2,3,8,10,12,18,17)]
#colnames(tab2) <- c("SOMAmer", "Chr","TSS","SNP","SNP_Pos38","Rank","Beta","Pval")
#tmp <- allpQTL[,c(1,8,10)]; colnames(tmp) <- c("SOMAmer","TopSNP","TopSNP_Pos38")
#tab2 <- inner_join( tab2, tmp, by="SOMAmer")
#tab2 <- tab2[,c(1:3,9:10,4:8)]
#tab2 <- inner_join(tab2, annota[,c("seqid_in_sample","uniprot_id","target","entrezgenesymbol","targetfullname")],
#                   by=c("SOMAmer"="seqid_in_sample"))
#tab2 <- tab2[,c("SOMAmer", "target","targetfullname","uniprot_id","entrezgenesymbol","Chr","TSS","TopSNP","TopSNP_Pos38","SNP","SNP_Pos38","Rank","Beta","Pval")]
#write_tsv(tab2,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/tab2.txt")

tab2  <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/tab2.txt")

SOMA <- unique(tab2$SOMAmer)
library( snpStats)

R2 <- numeric()
  tmp <- tab2[tab2$SOMAmer==SOMA[i],]
  for (j in 1:nrow(tmp)){
    if(tmp$TopSNP[1]  == tmp$SNP[j]){
      R2 <- c(R2, 1)
    }else{
      a <- read.plink(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr",tmp$Chr[1]),select.snps=c(tmp$TopSNP[1], tmp$SNP[j]))
      a <- as(a$genotypes, Class="numeric")
      if(sum(is.na(a))!=0){
        a[,1][is.na(a[,1])] <- mean(a[,1], na.rm = T)
        a[,2][is.na(a[,2])] <- mean(a[,2], na.rm = T)
      }
      R2 <- c(R2, cor(a[,1],a[,2])^2)
    }
  }

saveRDS(R2, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/tab2/R2_", i, ".rds"))
