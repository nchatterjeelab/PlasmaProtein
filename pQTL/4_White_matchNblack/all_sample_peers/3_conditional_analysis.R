
rm(list=ls())


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/4_White_matchNblack/all_sample_peers/3_conditional_analysis"))


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/vcf_file/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/conditional"))

library(readr)
library(stringr)


a <- paste0("#!/usr/bin/env bash
#$ -N Wsubmit
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

")

for (r in 1:10){

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/vcf_file/",r))
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation/",r))
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/conditional/",r))


  b <- paste0("#!/usr/bin/env bash
#$ -N W",r,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load R/3.6.1

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \\
--threads 1 \\
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr$SGE_TASK_ID \\
--keep /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/White_matchNblack_ID_",r,".txt \\
--export vcf \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/vcf_file/",r,"/chr$SGE_TASK_ID

cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/vcf_file/",r,"
bgzip chr$SGE_TASK_ID.vcf && tabix -p vcf chr$SGE_TASK_ID.vcf.gz

cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/bed_file/",r,"/
bgzip chr$SGE_TASK_ID.bed && tabix -p bed chr$SGE_TASK_ID.bed.gz

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation/",r,"/chr$chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/vcf_file/",r,"/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/bed_file/",r,"/chr$SGE_TASK_ID.bed.gz \\
--permute 100 \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation/",r,"/chr$SGE_TASK_ID/permutations.txt

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation/",r,"/chr$SGE_TASK_ID/permutations.txt 0.05 \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation/",r,"/chr$SGE_TASK_ID/permutations_all

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/conditional/",r,"/chr$SGE_TASK_ID

/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/vcf_file/",r,"/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/bed_file/",r,"/chr$SGE_TASK_ID.bed.gz \\
--mapping /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/permutation/",r,"/chr$SGE_TASK_ID/permutations_all.thresholds.txt \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/conditional/",r,"/chr$SGE_TASK_ID/conditional.txt


")


  writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/4_White_matchNblack/all_sample_peers/3_conditional_analysis/', r, '.sh'))

}

