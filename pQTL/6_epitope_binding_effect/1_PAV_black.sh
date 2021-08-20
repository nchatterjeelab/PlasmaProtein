#!/usr/bin/env bash
#$ -N PAV_b
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

ethnic='Black'

runr(){
    R CMD BATCH --no-save --no-restore "$1"  1_PAV.R ${ethnic}.Rout
}
runr "--args ethnic='${ethnic}'"
