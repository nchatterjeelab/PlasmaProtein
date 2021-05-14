
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race='black'

res <- '/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/'

dir.create(paste0(res))
dir.create(paste0(res, 'peernum'))
dir.create(paste0(res, 'peernum/peerandresid'))
dir.create(paste0(res, 'peernum/other'))
dir.create(paste0(res, 'peernum/covariates'))
dir.create(paste0(res, 'peernum/invrankpheno'))


a <- character()
b <- character()
n <- 0
for (peer in 1:10){
  n <- n+1
  n_peer <- 10*peer
a <- paste0(a, "PEER_wb.",race,".",n_peer,".Rout
")
b <- paste0(b, "race='",race,"' n_peer='",n_peer,"'
")
}

writeLines(a,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/1_peer_calculation/ind_Rout.txt"))
writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/1_peer_calculation/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N PEER
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b <  ind_para.txt

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  PEER_wb.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")
writeLines(a, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/1_peer_calculation/PEER.sh"))

