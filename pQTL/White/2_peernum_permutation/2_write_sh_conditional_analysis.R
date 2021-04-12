
rm(list=ls())


dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/2_peernum_permutation')
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation'))

library(readr)
library(stringr)

for (n_peer in (0:25)*10){
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/2_peernum_permutation/',n_peer))
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/',n_peer))
for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N Wp",n_peer,"c", i, "
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load R/3.6.1

cd /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/bed_file
bgzip chr",i,".bed && tabix -p bed chr",i,".bed.gz

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/bed_file/chr",i,".bed.gz \\
--permute 100 \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations.txt

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations.txt 0.05 \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/bed_file/chr",i,".bed.gz \\
--mapping /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/",n_peer,"/conditional/chr",i,"/conditional.txt


")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/2_peernum_permutation/',n_peer,'/chr', i, '.sh'))

  
  print(paste0("p",n_peer,"c",i))
}

}

#### backup 1%

dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/2_peernum_permutation/backup0.01")
dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/")

for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N Wp120c", i, "
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load conda_R/4.0

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/permutation
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/permutation/chr",i,"

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/permutation/chr",i,"/permutations.txt 0.01 \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/permutation/chr",i,"/permutations_all

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/conditional
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/conditional/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/bed_file/chr",i,".bed.gz \\
--mapping /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/permutation/chr",i,"/permutations_all.thresholds.txt \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/backup0.01/conditional/chr",i,"/conditional.txt
")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/2_peernum_permutation/backup0.01/chr', i, '.sh'))
}



dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/2_peernum_permutation/backup0.01")
dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/")

for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N Wp70c", i, "
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load conda_R/4.0

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/permutation
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/permutation/chr",i,"

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/permutation/chr",i,"/permutations.txt 0.01 \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/permutation/chr",i,"/permutations_all

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/conditional
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/conditional/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/bed_file/chr",i,".bed.gz \\
--mapping /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/permutation/chr",i,"/permutations_all.thresholds.txt \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/70/backup0.01/conditional/chr",i,"/conditional.txt
")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/2_peernum_permutation/backup0.01/chr', i, '.sh'))
}