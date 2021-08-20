#!/usr/bin/env bash
#$ -N impute_prot_AA
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


module load R/3.6.1

R CMD BATCH --no-save --no-restore 1_1000G_imputed_protein_AA.R

