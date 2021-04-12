#!/usr/bin/env bash
#$ -N jointPWAS
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=10G
#$ -m e
#$ -t 1-1590
#$ -M jzhan218@jhu.edu

readarray -t prot < /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/protein.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1"  1_run.R ./1_run/${prot[$(($SGE_TASK_ID-1))]}
}
runr "--args prot='${prot[$(($SGE_TASK_ID-1))]}'"



