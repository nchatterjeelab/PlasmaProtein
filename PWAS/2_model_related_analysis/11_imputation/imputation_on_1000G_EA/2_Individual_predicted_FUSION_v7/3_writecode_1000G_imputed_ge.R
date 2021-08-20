################################################################
################################################################
################################################################

rm(list=ls())

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

a <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -m e
#$ -M jzhan218@jhu.edu

####################
")

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/11_imputation/imputation_on_1000G/2_Individual_predicted_FUSION_v7/codes")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v7")

for (tissue in tissue_list){

b <- paste0("#!/usr/bin/env bash
#$ -N ", tissue ,"
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../3_1000G_imputed_ge.R ",tissue,".Rout
}
runr \"--args tissue='",tissue,"' \"

")

  a <- paste0(a,
"
qsub ",tissue, ".sh
")


  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/11_imputation/imputation_on_1000G/2_Individual_predicted_FUSION_v7/codes/", tissue, ".sh"))

  print(tissue)
}

writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/11_imputation/imputation_on_1000G/2_Individual_predicted_FUSION_v7/codes/all.sh"))

