

rm(list=ls())


ethnic <- "White"

if(ethnic == "White"){
  n_peer <- 90
}else{
  n_peer <- 80
}


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/"))


b <- paste0("#!/usr/bin/env bash
#$ -N backup0.01_",ethnic,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -m e

module load R/3.6.1

mkdir /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/chr$SGE_TASK_ID

Rscript /dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/script/runFDR_cis_new.R \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum_permutation/",n_peer,"/permutation/chr$SGE_TASK_ID/permutations.txt 0.01 \\
/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/chr$SGE_TASK_ID/permutations_all


/dcl01/chatterj/data/jzhang2/TOOLS/QTLtools_1.2_CentOS7.8_x86_64/QTLtools_1.2_CentOS7.8_x86_64 cis \\
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/",ethnic,"/chr$SGE_TASK_ID.vcf.gz \\
--bed /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum_permutation/",n_peer,"/bed_file/chr$SGE_TASK_ID.bed.gz \\
--mapping /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/chr$SGE_TASK_ID/permutations_all.thresholds.txt \\
--window 500000 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/backup0.01/chr$SGE_TASK_ID/conditional.txt

")


writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/5_backup_fdr0.01/conditional_",ethnic,".sh"))





