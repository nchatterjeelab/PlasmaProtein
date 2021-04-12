
rm(list=ls())


dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_rerun')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_rerun')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp')

library(readr)
library(stringr)

seqid <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/seqid_autosomal_withSNP.txt")

library(stringr)
tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
tmp <- tmp[str_detect(tmp, "RDat")]
tmp <- gsub("[.].*", "", tmp)

seqid <- setdiff(seqid, tmp)


a <- "#!/usr/bin/env bash
#$ -N submitjobs_rerun
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

#SeqId_19584_33 # 1706
#SeqId_2979_8 # 1817
#SeqId_3623_84 # 1924

for(i in 1924:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -m e

module load R/3.6.1

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq/", seqid[i], " \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--save_hsq TRUE \\
--models enet


")
  
  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_rerun/', seqid[i], '.sh'))
  
  
a <- paste0(a,
"
qsub ../shfiles_rerun/", seqid[i],".sh
")
  
  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_rerun/submitjobs_rerun.sh')



########################


library(stringr)
tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs")
tmp <- tmp[str_detect(tmp, "e[0-9]")]
a <- character()
for (i in 1:length(tmp)){
  b <- readLines(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs/",tmp[i]))
  b <- paste(b,collapse="; ")
  if(length(b)==0)
    a[i] <- "NA"
  else
    a[i] <- b
  print(i)
}
which(str_detect(a,"core"))
tmp[which(str_detect(a,"core"))]
#SeqId_12627_97.e2151606


tmp1 <- tmp[which(str_detect(a,"core"))]
tmp_2 <-  gsub("[.].*", "", tmp1)
seqid <- tmp_2

#tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs")
#tmp <- tmp[str_detect(tmp, "e[0-9]")]
#a <- character()
#for (i in 1:length(tmp)){
#  b <- readLines(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs/",tmp[i]))
#  b <- paste(b,collapse="; ")
#  if(length(b)==0)
#    a[i] <- "NA"
#  else
#    a[i] <- b
#}
#which(str_detect(a,"core"))
#tmp[which(str_detect(a,"core"))]
#
#tmp1 <- tmp[which(str_detect(a,"core"))]
#tmp_1 <-  gsub("[.].*", "", tmp1)
#seqid <- tmp_2
##SeqId_18286_3.e1689015

a <- "#!/usr/bin/env bash
#$ -N submitjobs_rerun
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -m e

module load R/3.6.1

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq/", seqid[i], " \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--save_hsq TRUE \\
--models enet


")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_rerun/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles_rerun/", seqid[i],".sh
")

  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_rerun/submitjobs_rerun.sh')


################################################
################################################

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

for job in {2151312..2151561}
do
qdel ${job}
done

