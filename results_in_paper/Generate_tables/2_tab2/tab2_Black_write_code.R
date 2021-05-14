
################################################################
################################################################
################################################################

rm(list=ls())

ethnic <- "Black"

dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/tab2/"))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/1_Generate_tables/2_tab2/",ethnic))

tab2  <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/tab2.txt"))

SOMA <- unique(tab2$SOMAmer)
length(SOMA)

a <- paste0("#!/usr/bin/env bash
#$ -N R2_", ethnic,"
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
")

for (i in 1:length(SOMA)){


b <- paste0("#!/usr/bin/env bash
#$ -N ", ethnic,"_", i,"
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=5G
#$ -m e

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../tab2_",ethnic,".R R2_",i,".Rout
}
runr \"--args i='",i,"'\"

")


a <- paste0(a,
"
qsub R2_", i,".sh
")


    print(i)
    writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/1_Generate_tables/2_tab2/",ethnic, "/R2_", i,".sh"))

}
writeLines(a,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/1_Generate_tables/2_tab2/",ethnic,"/all.sh"))

