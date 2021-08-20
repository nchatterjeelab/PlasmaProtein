
rm(list=ls())

ethnic="White"

library(readr)
library(stringr)

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/1_training/")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/1_training/shfiles")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/1_training/submitjobs")
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/PWAS"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/PWAS/coefs"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/PWAS/tmp"))
"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/1_training/submitjobs/"

seqid <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/seqid_autosomal_withSNP.txt"))

if(ethnic == "White"){
  n_peer <- 90
}else{
  n_peer <- 80
}

a <- paste0("#!/usr/bin/env bash
#$ -N ", ethnic,"_PWAS
#$ -cwd
#$ -t 1-",length(seqid),"
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e

module load old_conda_R/3.6

readarray -t a < /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/seqid_autosomal_withSNP.txt

Rscript /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights_plinkthreads.R \\
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--PATH_gcta /dcs04/nilanjan/data/jzhang2/TOOLS/gcta_1.93.2beta/gcta64 \\
--bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/byseq_remove_ambiguous_snp/${a[$(($SGE_TASK_ID-1))]} \\
--verbose 2 \\
--models enet \\
--save_hsq TRUE \\
--tmp /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/PWAS/tmp/${a[$(($SGE_TASK_ID-1))]} \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/peernum/invrankpheno/",n_peer,"/${a[$(($SGE_TASK_ID-1))]}.pheno \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/PWAS/coefs/${a[$(($SGE_TASK_ID-1))]}
")
# --hsq_set 1 \\
# --models top1,lasso,enet \\

writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/1_training/submitjobs/ALL_",ethnic,".sh"))
