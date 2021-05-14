#!/bin/bash
#SBATCH --job-name=imp_prot
#SBATCH --time=72:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

module load R/3.6.1

R CMD BATCH --no-save --no-restore 1_1000G_imputed_protein.R

