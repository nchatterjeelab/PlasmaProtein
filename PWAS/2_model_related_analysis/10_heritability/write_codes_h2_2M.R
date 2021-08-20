
rm(list=ls())

ethnic <- "White"

if(ethnic == "White"){
  n_peer <- 90
}else{
  n_peer <- 80
}


dir.create('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/10_heritability')
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/10_heritability/',ethnic))

dir.create('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/h2_2M')
dir.create('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/tmp_2M')

seqid <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/seqid_autosomal_withSNP.txt"))


b <- paste0("#!/usr/bin/env bash
#$ -N ", ethnic,"_h2
#$ -cwd
#$ -t 1-",length(seqid),"
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e

module load old_conda_R/3.6

readarray -t a < /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/seqid_autosomal_withSNP.txt

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/compute_h2.R \\
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--PATH_gcta /dcs04/nilanjan/data/jzhang2/TOOLS/gcta_1.93.2beta/gcta64 \\
--bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/window2M/byseq_remove_ambiguous_snp/${a[$(($SGE_TASK_ID-1))]} \\
--verbose 2 \\
--tmp /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/tmp_2M/${a[$(($SGE_TASK_ID-1))]} \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/h2_2M/${a[$(($SGE_TASK_ID-1))]} \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum/invrankpheno/", n_peer, "/${a[$(($SGE_TASK_ID-1))]}.pheno \\
--save_hsq TRUE \\
--models enet

")


writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/10_heritability/',ethnic,'/submit.sh'))

