
rm(list=ls())


dir.create('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/7_AASK/2_peernum_permutation/submit')
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK'))
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation'))

library(readr)
library(stringr)


a <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -m e
")

for (n_peer in (0:10)*10){
  dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/',n_peer))

b <- paste0("#!/usr/bin/env bash
#$ -N Bp",n_peer,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=20G,h_vmem=20G,h_fsize=30G
#$ -m e

module load R/3.6.1

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation
mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/vcf/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--permute 100 \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
 /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt 0.05 \\
 /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations_all

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional
mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/vcf/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--mapping /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations_all.thresholds.txt \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/conditional/chr$SGE_TASK_ID/conditional.txt


")

  writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/7_AASK/2_peernum_permutation/submit/',n_peer,'.sh'))

  a <- paste0(a, "qsub ",n_peer,".sh
")

  print(paste0("p",n_peer))

}

writeLines(a,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/7_AASK/2_peernum_permutation/submit/ALL.sh'))

