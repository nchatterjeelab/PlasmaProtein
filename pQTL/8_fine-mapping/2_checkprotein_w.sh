#!/usr/bin/env bash
#$ -N wLD
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

R CMD BATCH --no-save --no-restore 2_checkprotein_w.R

