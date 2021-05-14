#######################################
# PEER

rm(list=ls())

library(readr)

args <- commandArgs(T)

for(i in 1:length(args)){
     eval(parse(text=args[[i]]))
}

n_peer <- as.integer(n_peer)

soma <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/prot.txt")

res <- '/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/'

################

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt')

#################################################################

## covariates

library(dplyr)

covariates <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/aask.cov')
write_tsv(covariates, paste0(res, 'peernum/covariates/covariates.', n_peer, '.cov'))

#################################################################

# inv-rank phenotype

library(RNOmni)

dir.create(paste0(res, 'peernum/invrankpheno/', n_peer))
dir.create(paste0(res, 'peernum/results/'))
dir.create(paste0(res, 'peernum/results/', n_peer))

for(i in 1:length(seqid)){

  dat <- data.frame(y=soma[[seqid[i]]], covariates[,-1:-2])
  fit <- lm(y~., dat)
  dat <- data.frame(FID=0,IID=soma$gene_id, y=rankNorm(residuals(fit))); colnames(dat)[3] <- seqid[i]
  write_tsv(dat, col_names=F, paste0(res, 'peernum/invrankpheno/', n_peer, '/', seqid[i], '.pheno'))

  print(i)

}


