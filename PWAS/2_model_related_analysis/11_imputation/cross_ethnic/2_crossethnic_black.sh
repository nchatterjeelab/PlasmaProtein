#!/usr/bin/env bash
#$ -N cross_b
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=10G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

R CMD BATCH --no-save --no-restore 2_crossethnic_black.R





