
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(dplyr)
library(stringr)
library(genio)

AA <- read_tsv("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_AA_hg38.pos")
EA <- read_tsv("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_EA_hg38.pos")
protein1 <- intersect(AA$WGT,EA$WGT) # 1096
protein2 <- union(AA$WGT,EA$WGT)
protein <- unique(c(protein1,protein2))
protein <- gsub(".wgt.RDat","",protein)
writeLines(protein, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/protein.txt")



a <- paste0("#!/usr/bin/env bash
#$ -N LD
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -m e
#$ -t 1-",length(protein),"
#$ -M jzhan218@jhu.edu

readarray -t prot < /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/protein.txt

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--threads 1 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq_remove_ambiguous_snp/${prot[$(($SGE_TASK_ID-1))]} \\
--allow-no-sex \\
--r bin4 \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/LD/White/${prot[$(($SGE_TASK_ID-1))]}

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--threads 1 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/${prot[$(($SGE_TASK_ID-1))]} \\
--allow-no-sex \\
--r bin4 \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/LD/Black/${prot[$(($SGE_TASK_ID-1))]}

")
writeLines(a, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/jointPWAS/LD.sh"))

