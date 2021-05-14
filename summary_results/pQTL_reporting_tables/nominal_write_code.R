
rm(list=ls())

library(readr)
library(stringr)

ethnic <- "Black"

if(ethnic == "White"){
  n_peer <- 120
}else{
  n_peer <- 70
}

dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Generate_reporting_tables/nominal"))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/nominal"))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Generate_reporting_tables/nominal/",ethnic))

for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N nomi_", i, "
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=30G,chatterjee
#$ -m e

module load R/3.6.1

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/",ethnic,"/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum_permutation/",n_peer,"/bed_file/chr",i,".bed.gz \\
--nominal 1 \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Generate_reporting_tables/nominal/",ethnic,"/chr",i,".txt

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Generate_reporting_tables/nominal/',ethnic,'/chr', i, '.sh'))

}

ethnic <- "Black"
for (i in 1:22){
  system(paste0("mv /dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Generate_reporting_tables/nominal/",ethnic,"/chr",i,".txt /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/nominal/"))
  print(i)
}
