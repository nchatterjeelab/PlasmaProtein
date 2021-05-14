################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(stringr)

tissue_list <- readLines("/home-2/jzhan218@jhu.edu/work/jzhan218/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

a <- paste0("#!/bin/bash
#SBATCH --job-name=submit
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

####################
")

dir.create("/home-2/jzhan218@jhu.edu/work/jzhan218/PWAS/codes/White/1_invrank/3_PWAS/3_conditonal_analysis_PrediXcan/2_Individual_predicted_FUSION/codes/")
for (tissue in tissue_list){

b <- paste0("#!/bin/bash
#SBATCH --job-name=", tissue ,"
#SBATCH --time=72:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../3_1000G_imputed_ge.R ",tissue,".Rout
}
runr \"--args tissue='",tissue,"' \"


")

  a <- paste0(a,
"
sbatch ",tissue, ".sh
")


  writeLines(b,  paste0("/home-2/jzhan218@jhu.edu/work/jzhan218/PWAS/codes/White/1_invrank/3_PWAS/3_conditonal_analysis/2_Individual_predicted_FUSION/codes/", tissue, ".sh"))

  print(tissue)
}

writeLines(a,  paste0("/home-2/jzhan218@jhu.edu/work/jzhan218/PWAS/codes/White/1_invrank/3_PWAS/3_conditonal_analysis/2_Individual_predicted_FUSION/codes/all.sh"))


