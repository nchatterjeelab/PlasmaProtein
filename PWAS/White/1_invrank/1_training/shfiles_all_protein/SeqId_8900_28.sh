#!/usr/bin/env bash
#$ -N SeqId_8900_28
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=10G
#$ -m e

module load R/3.6.1

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights.R \
--PATH_plink /users/jzhang2/RESEARCH/tools/plink/plink2 \
--PATH_gcta /users/jzhang2/RESEARCH/tools/gcta_1.93.0beta/gcta64 \
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq/SeqId_8900_28 \
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_8900_28 \
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/coefs/SeqId_8900_28 \
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/90/SeqId_8900_28.pheno \
--models enet



/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq/SeqId_8900_28.fam