
rm(list=ls())


dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF")


for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N chr", i, "
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -pe local 3
#$ -m e

module load R/3.6.1

/users/jzhang2/RESEARCH/tools/plink/plink2 \\
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr",i," \\
--freq \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/MAF/chr",i,"

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/8_MAF/chr', i, '.sh'))


  print(i)
}

