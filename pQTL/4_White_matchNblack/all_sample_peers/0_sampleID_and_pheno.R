#######################################



rm(list=ls())

library(readr)
library(dplyr)
library(stringr)
library(RNOmni)

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')

set.seed(1)

dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/'))
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/'))
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/invrankpheno'))

for (r in 1:10){

  print(paste0("r=",r))

  dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/invrankpheno/',r))

  White_matchNblack <- sort(sample(1:7213, 1871))
  writeLines(as.character(White_matchNblack), paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/White_matchNblack_",r,".txt"))
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/90/SeqId_10000_28.pheno"), stringsAsFactors = F)
  random_sample <- tmp$V2[White_matchNblack]
  writeLines(random_sample, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/White_matchNblack_ID_",r,".txt"))

  covariates <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/covariates/covariates.90.cov'), col_types = cols())
  covariates <- covariates[covariates$IID %in% random_sample, ]

  for(i in 1:length(seqid)){
    gene <-  seqid[i]
    #gene <- "SeqId_6919_3"
    tmp <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pheno/', gene, '.pheno'), col_types = cols())
    dat <- inner_join(tmp,covariates,by = c("FID", "IID"))
    tmp <- dat[,1:2]; dat <- dat[,-1:-2]; colnames(dat)[1] <- "y"
    fit <- lm(y~., dat)
    dat <- data.frame(tmp, y=RankNorm(residuals(fit))); colnames(dat)[3] <- gene
    write_tsv(dat, col_names=F, paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/invrankpheno/',r,'/', gene, '.pheno'))

    print(i)
  }

}


#for r in {1..10}
#do
#mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/Random_10times/summary_data/${r}
#/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \
#--threads 1 \
#--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/byseq/SeqId_6919_3 \
#--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/Random_10times/invrankpheno/${r}/SeqId_6919_3.pheno \
#--glm allow-no-covars \
#--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/Random_10times/summary_data/${r}/SeqId_6919_3.pheno
#done



