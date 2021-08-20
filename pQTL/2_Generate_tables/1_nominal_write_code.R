
rm(list=ls())

library(readr)
library(stringr)

ethnic <- "Black"

if(ethnic == "White"){
  n_peer <- 90
}else{
  n_peer <- 80
}

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/2_Generate_tables/nominal"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/nominal"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/2_Generate_tables/nominal/",ethnic))


b <- paste0("#!/usr/bin/env bash
#$ -N nomi_",ethnic,"
#$ -t 1-22
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=30G,chatterjee
#$ -m e

module load R/3.6.1

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/",ethnic,"/chr$SGE_TASK_ID.vcf.gz \\
/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--nominal 1 \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/nominal/chr$SGE_TASK_ID.txt

")

writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/2_Generate_tables/nominal/',ethnic,'/nomi.sh'))

