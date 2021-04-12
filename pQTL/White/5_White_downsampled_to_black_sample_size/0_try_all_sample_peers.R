#######################################

library(readr)
library(dplyr)

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')

random_sample <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/White_matchNblack_ID.txt")

n_peer <- 120
covariates <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/covariates/covariates.', n_peer, '.cov'))
covariates <- covariates[covariates$IID %in% random_sample, ]

library(RNOmni)
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/invrankpheno/'))
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/invrankpheno/', n_peer))

for(i in 1:length(seqid)){

  tmp <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pheno/', seqid[i], '.pheno') )

  dat <- inner_join(tmp,covariates)
  tmp <- dat[,1:2]; dat <- dat[,-1:-2]; colnames(dat)[1] <- "y"

  fit <- lm(y~., dat)
  dat <- data.frame(tmp, y=rankNorm(residuals(fit))); colnames(dat)[3] <- seqid[i]
  write_tsv(dat, col_names=F, paste0(res, 'peernum/invrankpheno/', n_peer, '/', seqid[i], '.pheno'))

  print(i)

}


for(i in 1:length(seqid)){

  tmp <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pheno/', seqid[i], '.pheno') )
  tmp <- inner_join(tmp, covariates, by = c("FID", "IID"))
  dat <- tmp[,-1:-2]; colnames(dat)[1] <- "y"
  fit <- lm(y~., dat)
  dat <- data.frame(tmp[,1:2], y=rankNorm(residuals(fit))); colnames(dat)[3] <- seqid[i]
  write_tsv(dat, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/invrankpheno/try',r,'/', seqid[i], '.pheno'))

  print(i)

}


