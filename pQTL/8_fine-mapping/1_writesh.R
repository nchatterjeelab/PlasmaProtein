
################################################################
################################################################
################################################################

rm(list=ls())

ethnic='Black'
library(readr)

seqid <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/seqid_autosomal_withSNP.txt"))


if(ethnic == "White"){
  n_peer <- 90
}else{
  n_peer <- 80
}

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping"))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/8_fine-mapping"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/8_fine-mapping/submitjobs"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping/LD"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping/summary_stat"))


b <- paste0("#!/usr/bin/env bash
#$ -N ", ethnic,"
#$ -t 1-",length(seqid),"
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -cwd
#$ -m e

readarray -t a < /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/seqid_autosomal_withSNP.txt

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--threads 1 \\
--bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/byseq/${a[$(($SGE_TASK_ID-1))]} \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/peernum/invrankpheno/", n_peer, "/${a[$(($SGE_TASK_ID-1))]}.pheno \\
--glm cols=+a1freq,-a1countcc allow-no-covars --allow-no-sex \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping/summary_stat/${a[$(($SGE_TASK_ID-1))]}

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink \\
--threads 1 \\
--keep-allele-order \\
--bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/byseq/${a[$(($SGE_TASK_ID-1))]} \\
--r bin4 \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping/LD/${a[$(($SGE_TASK_ID-1))]}

")

writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/8_fine-mapping/submitjobs/ALL_",ethnic,".sh"))
paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/8_fine-mapping/submitjobs/ALL_",ethnic,".sh")


for(i in 1:length(seqid)){
  system(paste0("cp /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/window1M/byseq/", seqid[i],".bim /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/", ethnic, "/pQTL/fine-mapping/LD/"))
  print(i)
}
