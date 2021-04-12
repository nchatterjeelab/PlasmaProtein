#!/usr/bin/env bash
#$ -N v3w230
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

r=white
s='230'

runr(){
    R CMD BATCH --no-save --no-restore "$1" PEER_wb.R PEER_wb.${r}.${s}.Rout
}
runr "--args race='${r}' n_peer=${s}"
