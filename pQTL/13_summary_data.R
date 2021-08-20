### Summary data for cis-heritable genes

##################################################
## clean summary data
##################################################

rm(list=ls())

library(readr)

a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos")

seqid <- gsub(".wgt.RDat","",a$WGT)

annota <- read_tsv('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
annota <- annota[match(seqid, annota$seqid_in_sample),]

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat_website")

a <- paste0("#!/usr/bin/env bash
#$ -N clean_White
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat
####################
")

for(i in 1:length(seqid)){ a <- paste0(a, "
cp ",seqid[i],".PHENO1.glm.linear ../summary_stat_website/",annota$entrezgenesymbol[i],".", seqid[i], ".txt
") }
writeLines(a,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/13_summary_data/clean_White.sh'))


##################################################


rm(list=ls())

library(readr)

a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_AA_hg38.pos")

seqid <- gsub(".wgt.RDat","",a$WGT)

annota <- read_tsv('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
annota <- annota[match(seqid, annota$seqid_in_sample),]

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat_website")

a <- paste0("#!/usr/bin/env bash
#$ -N clean_Black
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
cd /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat
####################
")

for(i in 1:length(seqid)){ a <- paste0(a, "
cp ",seqid[i],".PHENO1.glm.linear ../summary_stat_website/",annota$entrezgenesymbol[i],".", seqid[i], ".txt
") }
writeLines(a,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/13_summary_data/clean_Black.sh'))



