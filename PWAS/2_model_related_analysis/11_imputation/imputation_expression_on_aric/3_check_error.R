################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(stringr)

tissue_list <- readLines("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

for (tissue in tissue_list){
  tmp <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/11_imputation/imputation_expression_on_aric/codes/",tissue,".Rout"))
  if(tmp[length(tmp)-2] != "> proc.time()"){
    print(tissue)
  }
}


a <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
")

for (tissue in c("Adipose_Subcutaneous",
"Brain_Amygdala",
"Brain_Cerebellum",
"Brain_Hypothalamus",
"Brain_Nucleus_accumbens_basal_ganglia",
"Testis")){

b <- paste0("#!/usr/bin/env bash
#$ -N ", tissue, "
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../1_imputed_ge_on_aric.R ",tissue,".Rout
}
runr \"--args tissue='",tissue,"' \"

")

  a <- paste0(a,
"
qsub ",tissue, ".sh
")


  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/11_imputation/imputation_expression_on_aric/codes/", tissue, ".sh"))

  print(tissue)
}

writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/11_imputation/imputation_expression_on_aric/codes/all_rerun.sh"))


