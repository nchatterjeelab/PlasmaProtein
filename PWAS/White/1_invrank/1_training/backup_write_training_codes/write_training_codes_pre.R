
rm(list=ls())


dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/outfiles')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp')

library(readr)
library(stringr)

seqid <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/seqid_autosomal_withSNP.txt")

a <- "#!/usr/bin/env bash
#$ -N submitjobs
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=20G
#$ -m e

module load R/3.6.1

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq/", seqid[i], " \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--save_hsq TRUE \\
--models enet


")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles/", seqid[i],".sh
")

  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs/submitjobs.sh')


a <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs") # coefs_1
b <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun_1") # coefs_2
c <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun") # coefs_3

length(a) #2389
length(b) #1219
length(c) #51

length(intersect(a,b)) #271
length(intersect(b,c)) #1
length(intersect(a,c)) #30

d <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs") # 3412
