########################################################################
## keep non-zero weights only

a <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs/")
a <- a[stringr::str_detect(a,"RDat")]

rec0 <- character()
for (i in 1:length(a)) {
    load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs/",a[i]))
    m <- wgt.matrix[,"enet"]!=0
    if(sum(m)==0){
        rec0 <- c(rec0, a[i])
        next
    }else{
       wgt.matrix <- wgt.matrix[m,,drop=FALSE]
       snps <- snps[m,,drop=FALSE]
       rownames(snps) <- NULL
    }
    
    save(cv.performance, hsq, hsq.pv, N.tot, snps, wgt.matrix,
        file=paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA/",a[i]))

    print(i)

}


a <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/coefs/")
a <- a[stringr::str_detect(a,"RDat")]

rec0 <- character()
for (i in 1:length(a)) {
    load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/coefs/",a[i]))
    m <- wgt.matrix[,"enet"]!=0
    if(sum(m)==0){
        rec0 <- c(rec0, a[i])
        next
    }else{
       wgt.matrix <- wgt.matrix[m,,drop=FALSE]
       snps <- snps[m,,drop=FALSE]
       rownames(snps) <- NULL
    }
    
    save(cv.performance, hsq, hsq.pv, N.tot, snps, wgt.matrix,
        file=paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_AA/",a[i]))

    print(i)

}


rec0 ## EA
# "SeqId_11140_56.wgt.RDat" "SeqId_12963_1.wgt.RDat" 


rec0 ## AA
# [1] "SeqId_10073_22.wgt.RDat" "SeqId_11289_31.wgt.RDat"
# [3] "SeqId_14131_37.wgt.RDat" "SeqId_16831_7.wgt.RDat" 
# [5] "SeqId_17685_9.wgt.RDat"  "SeqId_17750_8.wgt.RDat" 
# [7] "SeqId_4673_13.wgt.RDat"  "SeqId_4931_59.wgt.RDat" 
# [9] "SeqId_5070_76.wgt.RDat" 

rec0 <- c("SeqId_11140_56.wgt.RDat","SeqId_12963_1.wgt.RDat" )
a <- readr::read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos"))
a <- a[!(a$WGT %in% rec0),]
readr::write_tsv(a, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos"))

rec0 <- c("SeqId_10073_22.wgt.RDat","SeqId_11289_31.wgt.RDat",
"SeqId_14131_37.wgt.RDat","SeqId_16831_7.wgt.RDat",
"SeqId_17685_9.wgt.RDat","SeqId_17750_8.wgt.RDat",
"SeqId_4673_13.wgt.RDat","SeqId_4931_59.wgt.RDat",
"SeqId_5070_76.wgt.RDat")
a <- readr::read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_AA_hg38.pos"))
a <- a[!(a$WGT %in% rec0),]
readr::write_tsv(a, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_AA_hg38.pos"))


#############################################
### PWAS

library(dplyr)
library(readr)
results <- tibble()
for (chr in 1:22) {
  results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/PWAS_CI/chr", chr, ".out")))
  if(chr==6){
    results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/PWAS_CI/chr6.out.MHC")))
  }
}
write_tsv(results, "/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/PWAS_CI.out")


########################################################################

# TWAS

rm(list=ls())
dir.create("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/submit/Gout/TWAS_CI")

tissue_list <- readLines("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/GTEx_V7_tissue_list.txt")

b <- paste0("#!/usr/bin/env bash
#$ -N Gout
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (tissue in tissue_list){
dir.create(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI/",tissue))

a <- paste0("#!/usr/bin/env bash
#$ -N Gout_",tissue,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript /dcl01/chatterj/data/jzhang2/PWAS_tutorial/scripts/PWAS.assoc_test_CI.R \\
--sumstats /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Gout_EA_cleaned_rsid_BETA_SE.txt \\
--weights /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/",tissue,".P01.pos \\
--weights_dir /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/ \\
--ref_ld_chr /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/LDREF/1000G.EUR. \\
--force_model enet \\
--chr $SGE_TASK_ID \\
--out /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI/",tissue,"/chr$SGE_TASK_ID.out

")
writeLines(a, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/submit/Gout/TWAS_CI/",tissue,".sh"))

  b <- paste0(b, "
qsub ",tissue,".sh")

}
writeLines(b, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/submit/Gout/TWAS_CI/ALL.sh"))




