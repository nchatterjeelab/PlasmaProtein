
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
ethnic <- "Black"

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/tab2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/2_Generate_tables/2_tab2/",ethnic))

tab2  <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/tab2.txt"))

SOMA <- unique(tab2$SOMAmer)
length(SOMA)

b <- paste0("#!/usr/bin/env bash
#$ -N ", ethnic,"
#$ -cwd
#$ -t 1-",length(SOMA),"
#$ -l mem_free=2G,h_vmem=2G,h_fsize=5G
#$ -m e

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../tab2_",ethnic,".R R2_$SGE_TASK_ID.Rout
}
runr \"--args i='$SGE_TASK_ID'\"

")

writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/2_Generate_tables/2_tab2/",ethnic, "/R2.sh"))


