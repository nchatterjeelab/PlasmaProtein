#######################################
# PEER

rm(list=ls())

n_peer <- 0

res <- '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/'

dir.create(paste0(res))
dir.create(paste0(res, 'peernum'))
dir.create(paste0(res, 'peernum/peerandresid'))
dir.create(paste0(res, 'peernum/other'))
dir.create(paste0(res, 'peernum/covariates'))
dir.create(paste0(res, 'peernum/invrankpheno'))
dir.create(paste0(res, 'peernum/invrankpheno/', n_peer))
dir.create(paste0(res, 'peernum/results/'))
dir.create(paste0(res, 'peernum/results/', n_peer))


## covariates
library(readr)
library(RNOmni)
library(fastDummies)

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')

aric.cov <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/aric.cov', header=TRUE, stringsAsFactors=F)
pc <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pc.txt', header=TRUE, stringsAsFactors=F)
pc <- pc[,-1]

covariates <- cbind(aric.cov, pc)
aric.center <- aric.cov$v3center; aric.center <- model.matrix(~aric.center)
covariates <- cbind(covariates[1:4], aric.center[,-1], covariates[,6:ncol(covariates)])

write_tsv(covariates, paste0(res, 'peernum/covariates/covariates.', n_peer, '.cov'))


for(i in 1:length(seqid)){

  tmp <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pheno/', seqid[i], '.pheno') )
  dat <- data.frame(y=tmp[[3]], covariates[,-1:-2])
  fit <- lm(y~., dat)
  dat <- data.frame(tmp[,1:2], y=rankNorm(residuals(fit))); colnames(dat)[3] <- seqid[i]
  write_tsv(dat, paste0(res, 'peernum/invrankpheno/', n_peer, '/', seqid[i], '.pheno'))

  print(i)
}


