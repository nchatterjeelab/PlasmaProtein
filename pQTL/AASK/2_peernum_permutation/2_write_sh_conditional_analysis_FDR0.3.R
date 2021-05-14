
rm(list=ls())


dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/2_peernum_permutation_FDR0.3')
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation'))

library(readr)
library(stringr)

  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/'))


a <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -m e
")

for (n_peer in (0:10)*10){
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/2_peernum_permutation_FDR0.3/',n_peer))
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/',n_peer))
for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N Bp",n_peer,"c", i, "
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=30G
#$ -m e

module load R/3.6.1

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/permutation
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/permutation/chr",i,"

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
 /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations.txt 0.3 \\
 /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/permutation/chr",i,"/permutations_all

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/conditional
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/conditional/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/vcf/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation/",n_peer,"/bed_file/chr",i,".bed.gz \\
--mapping /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum_permutation_FDR0.3/",n_peer,"/conditional/chr",i,"/conditional.txt


")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/2_peernum_permutation_FDR0.3/',n_peer,'/chr', i, '.sh'))

  a <- paste0(a,"
cd /dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/2_peernum_permutation_FDR0.3/",n_peer,"
qsub chr", i, ".sh")
  
  print(paste0("p",n_peer,"c",i))
}

}

writeLines(a,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/AASK/2_peernum_permutation_FDR0.3/submit_cond_analysis.sh'))

