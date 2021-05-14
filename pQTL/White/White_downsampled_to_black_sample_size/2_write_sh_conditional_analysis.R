
rm(list=ls())


dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/')
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack'))

library(readr)
library(stringr)

a <- paste0("#!/usr/bin/env bash
#$ -N Wsubmit
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

")

#for (n_peer in (0:12)*10){

n_peer <- 120

  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/',n_peer))
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/',n_peer))
for(i in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N Wp",n_peer,"c", i, "
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load R/3.6.1

/users/jzhang2/RESEARCH/tools/plink/plink2 \\
--threads 1 \\
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr",i," \\
--keep /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/White_matchNblack_ID.txt \\
--export vcf \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/vcf_file/chr",i,"

cd c",n_peer,"/vcf_file
bgzip chr",i,".vcf && tabix -p vcf chr",i,".vcf.gz

cd /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/bed_file
bgzip chr",i,".bed && tabix -p bed chr",i,".bed.gz

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/permutation
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/permutation/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/vcf_file/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/bed_file/chr",i,".bed.gz \\
--permute 100 \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/permutation/chr",i,"/permutations.txt

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/permutation/chr",i,"/permutations.txt 0.05 \\
/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/permutation/chr",i,"/permutations_all

mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/conditional
mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/conditional/chr",i,"

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 \\
cis --vcf /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/vcf_file/chr",i,".vcf.gz \\
--bed /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/bed_file/chr",i,".bed.gz \\
--mapping /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/permutation/chr",i,"/permutations_all.thresholds.txt \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/conditional/chr",i,"/conditional.txt

")

  a <- paste0(a, 'cd /dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/',n_peer,'
qsub chr', i, '.sh

')

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/',n_peer,'/chr', i, '.sh'))

  
  print(paste0("p",n_peer,"c",i))
}

#}

writeLines(a,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/all.sh'))

