#!/usr/bin/env bash
#$ -N replication
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

R CMD BATCH --no-save --no-restore replication.R

