#!/usr/bin/env bash
#$ -N prediction
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

R CMD BATCH --no-save --no-restore prediction.R

