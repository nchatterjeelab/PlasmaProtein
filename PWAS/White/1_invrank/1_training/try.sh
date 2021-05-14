#!/usr/bin/env bash
#$ -N SeqId_8900_28
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=10G
#$ -m e

module load R/3.6.1

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights.R \
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/SeqId_14151_4 \
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_14151_4 \
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/coefs_remove_ambiguous_snp/SeqId_14151_4 \
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/70/SeqId_14151_4.pheno \
--save_hsq TRUE \
--models enet

opt$PATH_plink <- "/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2"
opt$PATH_gcta <- "/dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64"
opt$bfile <- "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/SeqId_14151_4"
opt$tmp<- "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_14151_4"
opt$out<- "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/coefs_remove_ambiguous_snp/SeqId_14151_4"
opt$pheno<- "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/70/SeqId_14151_4.pheno"
opt$save_hsq<- TRUE
opt$models <- "enet"

### Estimating heritability
sh: line 1: 32093 Aborted


/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 --allow-no-sex \
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_14151_4 \
--make-grm-bin \
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/tmp/SeqId_14151_4




