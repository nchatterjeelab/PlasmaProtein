
################################################################
################################################################
################################################################

rm(list=ls())


for (chr in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N G_chr", chr,"
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcl01/chatterj/data/jzhang2/1000G/GRCh38/original_files

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--vcf ALL.chr",chr,"_GRCh38.genotypes.20170504.vcf.gz \\
--vcf-half-call m \\
--max-alleles 2 \\
--make-bed \\
--out chr",chr,"


")
  print(chr)

  writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/1_reformat1000G/chr", chr,".sh"))

}


for (chr in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N G_chr", chr,"
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcl01/chatterj/data/jzhang2/1000G/GRCh38/original_files

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--bfile ./transform_bfile/chr",chr," \\
--rm-dup exclude-mismatch \\
--make-bed \\
--out chr",chr,"

")
  print(chr)

  writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/1_reformat1000G/chr", chr,".sh"))

}


