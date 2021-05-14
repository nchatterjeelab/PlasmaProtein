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

prot <- soma[,seqid]

###############
# PEER

library(peer)

model=PEER()
PEER_setPhenoMean(model, as.matrix(prot))

PEER_setNk(model, n_peer)

PEER_setAdd_mean(model, TRUE)

PEER_update(model)
factors = PEER_getX(model)

colnames(factors) <- c('MeanFactor',paste('peer', 1:n_peer, sep=''))
peer <- data.frame(gene_id = soma$gene_id, factors[,-1])
write_tsv(peer, paste0(res, 'peernum/peerandresid/peers.', race, '_', n_peer, '.txt'))

residuals = PEER_getResiduals(model); colnames(residuals) <- colnames(prot)
resid <- data.frame(gene_id = soma$gene_id, residuals)
write_tsv(resid, paste0(res, 'peernum/peerandresid/resid.', race, '_', n_peer, '.txt'))

precision = PEER_getAlpha(model)
write.table(precision, paste0(res, 'peernum/other/precision.', n_peer, '.txt'), col.names = F, row.names = F)
weights = PEER_getW(model)
write.table(weights, paste0(res, 'peernum/other/weights.', n_peer, '.txt'), col.names = F, row.names = F)

#################################################################

## covariates

library(dplyr)

covariates <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/aask.cov')
peer <- read_tsv(paste0(res, 'peernum/peerandresid/peers.', race, '_', n_peer, '.txt'))
covariates <- cbind(covariates, peer[,-1])
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


