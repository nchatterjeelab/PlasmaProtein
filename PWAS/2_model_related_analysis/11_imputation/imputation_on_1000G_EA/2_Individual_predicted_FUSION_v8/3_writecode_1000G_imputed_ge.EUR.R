################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(stringr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")

a <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -m e
#$ -M jzhan218@jhu.edu

####################
")

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/6_conditonal_analysis/2_Individual_predicted_FUSION_v8/codes.EUR")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v8.EUR")

for (tissue in tissue_list){

b <- paste0("#!/usr/bin/env bash
#$ -N ", tissue ,"
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../3_1000G_imputed_ge.EUR.R ",tissue,".Rout
}
runr \"--args tissue='",tissue,"' \"

")

  a <- paste0(a,
"
qsub ",tissue, ".sh
")


  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/6_conditonal_analysis/2_Individual_predicted_FUSION_v8/codes.EUR/", tissue, ".sh"))

  print(tissue)
}

writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/6_conditonal_analysis/2_Individual_predicted_FUSION_v8/codes.EUR/all.sh"))




################################################################
################################################################
################################################################

## check error


library(stringr)
slurm <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/6_conditonal_analysis/2_Individual_predicted_FUSION_v8/codes.EUR/")
slurm <- slurm[str_detect(slurm,".Rout")]
for (i in 1:length(slurm)){
  a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/6_conditonal_analysis/2_Individual_predicted_FUSION_v8/codes.EUR/",slurm[i]))
  tmp <- sum(str_detect(a," No space left on device"))
  if(tmp>0){
    tmp <- gsub("https://gtexv8fusion.s3-us-west-1.amazonaws.com/ALL/GTEXv8.ALL.","",str_split(a[1]," ")[[1]][4])
    tmp <- gsub(".tar.gz","",tmp)
    print(tmp)
  }
}




#for job in {45301430..45304591}
#do
#  scontrol hold ${job}
#done
#
#for job in {45301430..45304591}
#do
#  scontrol release ${job}
#done
#
#for job in {45301430..45304591}
#do
#  scancel ${job}
#done
#
#for job in {1744820..1746747}
#do
#  qdel ${job}
#  done
