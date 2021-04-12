
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt')

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
annota <- annota[match(seqid, annota$seqid_in_sample),]

chr <- annota$chromosome_name
TSS <- annota$transcription_start_site


n_peer=70

dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/results/', n_peer))

dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/2_peernum/', n_peer))
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/2_peernum/', n_peer, '/shfiles'))
dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/2_peernum/', n_peer, '/submitjobs'))

a <- paste0("#!/usr/bin/env bash
#$ -N sub_", n_peer, "
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu


####################
")

for(i in 1:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -m e


#/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
#--threads 1 \\
#--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq/", seqid[i], " \\
#--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/", n_peer, "/", seqid[i], ".pheno \\
#--glm allow-no-covars \\
#--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/", seqid[i], "
#
#/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink \\
#--threads 1 \\
#--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq/", seqid[i], " \\
#--r bin4 \\
#--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/LD/", seqid[i], "

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--freq \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/", seqid[i], "

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/7_fine-mapping/shfiles/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles/", seqid[i],".sh
")

  print(paste0(n_peer, "_", i))
}

writeLines(a,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/Black/7_fine-mapping/submitjobs/submitjobs.sh'))




for(i in 1:length(seqid)){
  system(paste0("cp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq/", seqid[i],".bim /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/"))
  print(i)
}

# /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/SeqId_5620_13.PHENO1.glm.linear


# pam SeqId_5620_13