#!/usr/bin/env bash
#$ -N aaskgeno
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink/chr$SGE_TASK_ID \
--threads 1 \
--export vcf \
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/vcf/chr$SGE_TASK_ID

cd /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/vcf/

bgzip -c chr$SGE_TASK_ID.vcf > chr$SGE_TASK_ID.vcf.gz
tabix -p vcf chr$SGE_TASK_ID.vcf.gz
