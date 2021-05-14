
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)

for (chr in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N G_chr", chr,"
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcl01/chatterj/data/jzhang2/1000G/GRCh38/original_files

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr",chr,"_GRCh38.genotypes.20170504.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr",chr,"_GRCh38.genotypes.20170504.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr",chr,"_GRCh38_sites.20170504.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr",chr,"_GRCh38_sites.20170504.vcf.gz.tbi

")
  print(chr)

  writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/0_download1000G/chr", chr,".sh"))

}

