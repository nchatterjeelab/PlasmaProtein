

rm(list=ls())

dir.create('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/1_peernum_permutation')
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL'))
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation'))

library(readr)
library(stringr)

a <- paste0("#!/usr/bin/env bash
#$ -N Wp_submit
#$ -cwd
#$ -m e
")

for (n_peer in (0:25)*10){
  dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/',n_peer))

b <- paste0("#!/usr/bin/env bash
#$ -N Wp",n_peer,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load R/3.6.1

cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/bed_file
bgzip chr$SGE_TASK_ID.bed && tabix -p bed chr$SGE_TASK_ID.bed.gz

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation
mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--permute 100 \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt 0.05 \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations_all

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional
mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--mapping /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations_all.thresholds.txt \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr$SGE_TASK_ID/conditional.txt


")

  writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/1_peernum_permutation/',n_peer,'.sh'))

    a <- paste0(a, "qsub ",n_peer,".sh
")
  print(paste0("p",n_peer))

}

writeLines(a,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/1_peernum_permutation/ALL.sh'))






rm(list=ls())

dir.create('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/1_peernum_permutation')
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL'))
dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation'))

library(readr)
library(stringr)

paste0("cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/bed_file
bgzip chr",i,".bed && tabix -p bed chr",i,".bed.gz")

a <- paste0("#!/usr/bin/env bash
#$ -N Bp_submit
#$ -cwd
#$ -m e
")

for (n_peer in (0:20)*10){
  dir.create(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/',n_peer))

b <- paste0("#!/usr/bin/env bash
#$ -N Bp",n_peer,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load R/3.6.1

cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/bed_file
bgzip chr$SGE_TASK_ID.bed && tabix -p bed chr$SGE_TASK_ID.bed.gz

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation
mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--permute 100 \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt 0.05 \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations_all

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional
mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--mapping /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations_all.thresholds.txt \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/chr$SGE_TASK_ID/conditional.txt


")

  writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/1_peernum_permutation/',n_peer,'.sh'))

    a <- paste0(a, "qsub ",n_peer,".sh
")
  print(paste0("p",n_peer))

}

writeLines(a,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/1_peernum_permutation/ALL.sh'))


