#!/usr/bin/env bash
#$ -N geno_w
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr$SGE_TASK_ID \
--threads 1 \
--export vcf \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White/chr$SGE_TASK_ID

cd /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/White

bgzip -c chr$SGE_TASK_ID.vcf > chr$SGE_TASK_ID.vcf.gz
tabix -p vcf chr$SGE_TASK_ID.vcf.gz

