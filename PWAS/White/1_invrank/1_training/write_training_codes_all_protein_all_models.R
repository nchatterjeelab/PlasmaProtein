
rm(list=ls())

dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_all_protein_all_models')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_all_protein_all_models')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_all_protein_all_models')
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
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e

module load old_conda_R/3.6

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights_plinkthreads.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq_remove_ambiguous_snp/", seqid[i], " \\
--verbose 2 \\
--hsq_set 1 \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_all_protein_all_models/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--models top1,lasso,enet

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_all_protein_all_models/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles_all_protein_all_models/", seqid[i],".sh
")

  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_all_protein_all_models/submitjobs.sh')

