
rm(list=ls())

library(readr)
library(stringr)

ethnic <- "White"

dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/nominal"))
dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/nominal")


for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N nomi_", i, "
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=30G,chatterjee
#$ -m e

module load R/3.6.1

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/120/vcf_file/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/120/bed_file/chr",i,".bed.gz \\
--nominal 1 \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/nominal/chr",i,".txt

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/nominal/chr', i, '.sh'))

}
