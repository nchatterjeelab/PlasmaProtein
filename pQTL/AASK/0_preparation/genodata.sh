#!/usr/bin/env bash
#$ -N aaskgeno
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=1000G
#$ -pe local 5
#$ -m e
#$ -M jzhan218@jhu.edu

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcl02/leased/kidney/AASK/static/Genetics/05_imputed_Topmed/chr_$SGE_TASK_ID/chr$SGE_TASK_ID.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--threads 5 \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink/chr$SGE_TASK_ID

