
################################################################
################################################################
################################################################

rm(list=ls())

ethnic <- "EUR"
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/2_ethnics/",ethnic))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/1000G/GRCh38/",ethnic))

for (chr in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N ",ethnic,"_", chr,"
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=20G
#$ -m e

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--bfile /dcl01/chatterj/data/jzhang2/1000G/GRCh38/original_files/chr",chr," \\
--keep /dcl01/chatterj/data/jzhang2/1000G/1000G_",ethnic,"_ID.txt \\
--make-bed \\
--out /dcl01/chatterj/data/jzhang2/1000G/GRCh38/",ethnic,"/chr",chr,"

")
  print(chr)

  writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/2_1000Gethnics/",ethnic,"/chr",chr,".sh"))

}


