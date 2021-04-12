#######################################
# PEER

rm(list=ls())

library(readr)

args <- commandArgs(T)

for(i in 1:length(args)){
     eval(parse(text=args[[i]]))
}

n_peer <- as.integer(n_peer)

if(race=='white'){
  soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/soma_visit_3_log2_SMP.txt')
  race1 <- 'White'
  res <- '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/'
}

if(race=='black'){
  soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/soma_visit_3_log2_SMP.txt')
  race1 <- 'Black'
  res <- '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/'
}

dir.create(paste0(res))
dir.create(paste0(res, 'peernum'))
dir.create(paste0(res, 'peernum/peerandresid'))
dir.create(paste0(res, 'peernum/other'))
dir.create(paste0(res, 'peernum/covariates'))
dir.create(paste0(res, 'peernum/invrankpheno'))

################

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')

prot <- soma[,seqid]

# ###############
# # PEER
#
# library(peer)
#
# model=PEER()
# PEER_setPhenoMean(model, as.matrix(prot))
#
# PEER_setNk(model, n_peer)
#
# PEER_setAdd_mean(model, TRUE)
#
# PEER_update(model)
# factors = PEER_getX(model)
#
# colnames(factors) <- c('MeanFactor',paste('peer', 1:n_peer, sep=''))
# peer <- data.frame(SampleId = soma$SampleId, factors[,-1])
# write_tsv(peer, paste0(res, 'peernum/peerandresid/peers.', race, '_', n_peer, '.txt'))
#
# residuals = PEER_getResiduals(model); colnames(residuals) <- colnames(prot)
# resid <- data.frame(SampleId = soma$SampleId, residuals)
# write_tsv(resid, paste0(res, 'peernum/peerandresid/resid.', race, '_', n_peer, '.txt'))
#
# precision = PEER_getAlpha(model)
# write.table(precision, paste0(res, 'peernum/other/precision.', n_peer, '.txt'), col.names = F, row.names = F)
# weights = PEER_getW(model)
# write.table(weights, paste0(res, 'peernum/other/weights.', n_peer, '.txt'), col.names = F, row.names = F)

#################################################################

## covariates

library(dplyr)
# library(foreign)
# tmp <- data.frame(FID=0, IID=soma$SampleId, stringsAsFactors = F)
# aric.cov <- read.dta('/dcl01/chatterj/data/jzhang2/pwas/pipeline/ARIC/derive37.dta')
# aric.cov <- aric.cov[,c('id','v3age31','racegrp','gender','v3center')] # ?center
# aric.cov <- inner_join(tmp, aric.cov, by=c('IID'='id'))
# aric.cov$racegrp <- NULL
# aric.cov$gender <- ifelse(aric.cov$gender == 'F', 2, 1)
# write_tsv(aric.cov, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/', race1,'/aric.cov'))

covariates <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/', race1,'/aric.cov'))
pc <- read.table(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/',race1,'/pc.txt'), header=TRUE, stringsAsFactors = F)
covariates <- inner_join(covariates, pc, by=c("IID"="SampleId"))

peer <- read_tsv(paste0(res, 'peernum/peerandresid/peers.', race, '_', n_peer, '.txt'))
covariates <- cbind(covariates, peer[,-1])
aric.center <- covariates$v3center; aric.center <- model.matrix(~aric.center)
aric.center <- as_tibble(aric.center)
covariates <- cbind(covariates[1:4], as.matrix(aric.center[,-1]), covariates[,6:ncol(covariates)])

write_tsv(covariates, paste0(res, 'peernum/covariates/covariates.', n_peer, '.cov'))

#################################################################

# inv-rank phenotype

library(RNOmni)

dir.create(paste0(res, 'peernum/invrankpheno/', n_peer))
dir.create(paste0(res, 'peernum/results/'))
dir.create(paste0(res, 'peernum/results/', n_peer))

for(i in 1:length(seqid)){
  
  tmp <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/',race1,'/pheno/', seqid[i], '.pheno') )
  dat <- data.frame(y=tmp[[3]], covariates[,-1:-2])
  fit <- lm(y~., dat)
  dat <- data.frame(tmp[,1:2], y=rankNorm(residuals(fit))); colnames(dat)[3] <- seqid[i]
  write_tsv(dat, col_names=F, paste0(res, 'peernum/invrankpheno/', n_peer, '/', seqid[i], '.pheno'))
  
  print(i)

}


