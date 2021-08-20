

rm(list=ls())

library(readr)

"cp -r /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink /fastscratch/myscratch/jzhang2/AASK/"

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt')

annota <- annota[match(seqid,annota$seqid_in_sample),]

chr <- annota$chromosome_name
TSS <- annota$transcription_start_site

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/geno/window1M")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/geno/window1M/byseq")

b <- "#!/usr/bin/env bash
#$ -N window1M
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0(b,
"

/users/jzhang2/RESEARCH/tools/plink/plink2 \\
--threads 1 \\
--bfile /fastscratch/myscratch/jzhang2/AASK/plink/chr", chr[i], " \\
--chr ", chr[i], " --from-bp ", TSS[i]-500000," --to-bp ", TSS[i]+500000, " \\
--make-bed \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/geno/window1M/byseq/", seqid[i], "

")
 print(i)
}

writeLines(b,  '/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/geno/window1M/extractgene.sh')


