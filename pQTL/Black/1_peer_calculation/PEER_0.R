#######################################
# PEER

rm(list=ls())

library(readr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

n_peer <- as.integer(n_peer)

res <- '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/'

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

aric.cov <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/aric.cov', header=TRUE, stringsAsFactors=F)
pc <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pc.txt', header=TRUE, stringsAsFactors=F)
pc <- pc[,-1]

covariates <- cbind(aric.cov, pc)
aric.center <- aric.cov$v3center; aric.center <- model.matrix(~aric.center)
covariates <- cbind(covariates[1:4], aric.center[,-1], covariates[,6:ncol(covariates)])

write_tsv(covariates, paste0(res, 'peernum/covariates/covariates.', n_peer, '.cov'))


for(i in 1:length(seqid)){

  tmp <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pheno/', seqid[i], '.pheno') )
  dat <- data.frame(y=tmp[[3]], covariates[,-1:-2])
  fit <- lm(y~., dat)
  dat <- data.frame(tmp[,1:2], y=rankNorm(residuals(fit))); colnames(dat)[3] <- seqid[i]
  write_tsv(dat, paste0(res, 'peernum/invrankpheno/', n_peer, '/', seqid[i], '.pheno'))

  print(i)
}

#/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/0/SeqId_10001_7.pheno