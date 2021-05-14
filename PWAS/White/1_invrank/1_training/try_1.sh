#!/usr/bin/env bash
#$ -N try
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e
#$ -M jzhan218@jhu.edu

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 --allow-no-sex \
--threads 1 \
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_14151_4 \
--make-grm-bin \
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_14151_4


#--threads 1 \

