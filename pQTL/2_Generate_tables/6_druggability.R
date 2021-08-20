

## results from revision_500Kb/4_AASK

rm(list=ls())

library(readr)
library(bigreadr)
library(stringr)

pos0 <- read_tsv("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_EA_hg38.pos")$ID
pos0 <- c(pos0, read_tsv("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_AA_hg38.pos")$ID)
pos0 <- unique(pos0)

pos <- read_tsv("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/files_to_share/PWAS_tutorial/Plasma_Protein_EA_hg38.pos")$ID
pos <- c(pos, read_tsv("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/files_to_share/PWAS_tutorial/Plasma_Protein_AA_hg38.pos")$ID)
pos <- unique(pos)

writeLines(pos, "/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/files_to_share/PWAS_tutorial/all_genes_EA_and_AA.txt")

