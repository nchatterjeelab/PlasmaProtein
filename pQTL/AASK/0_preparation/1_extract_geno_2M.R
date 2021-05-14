
#rm(list=ls())
#
#library(readr)
#
#
#annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
#seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
#
#
#chr <- annota$chromosome_name
#TSS <- annota$transcription_start_site
#
#dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre")
#dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/byseq")
#
#b <- "#!/usr/bin/env bash
##$ -N window2M
##$ -cwd
##$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
##$ -m e
##$ -M jzhan218@jhu.edu
#
#####################
#"
#
#for(i in 1:length(seqid)){
#
#b <- paste0(b,
#"
#
#/users/jzhang2/RESEARCH/tools/plink/plink2 \\
#--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr", chr[i], " \\
#--chr ", chr[i], " --from-bp ", TSS[i]-1000000," --to-bp ", TSS[i]+1000000, " \\
#--make-bed \\
#--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/byseq/", seqid[i], "
#
#####################
#
#")
# print(i)
#}
#
#writeLines(b,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M_pre/extractgene.sh')

################################################################
################################################################
################################################################



rm(list=ls())

library(readr)

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt')

annota <- annota[match(seqid,annota$seqid_in_sample),]

chr <- annota$chromosome_name
TSS <- annota$transcription_start_site

dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/window2M")
dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/window2M/byseq")

b <- "#!/usr/bin/env bash
#$ -N window2M
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0(b,
"

/users/jzhang2/RESEARCH/tools/plink/plink2 \\
--threads 1 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink/chr", chr[i], " \\
--chr ", chr[i], " --from-bp ", TSS[i]-1000000," --to-bp ", TSS[i]+1000000, " \\
--make-bed \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/window2M/byseq/", seqid[i], "

")
 print(i)
}

writeLines(b,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/window2M/extractgene.sh')


########################################################################
########################################################################




############################################################################
############################################################################
### gene with snp info
#
#tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq_remove_ambiguous_snp")
#library(stringr)
#length(tmp[str_detect(tmp,"bed")]) #[1] 4660
#seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
#length(seqid) #[1] 4665
#tmp <- tmp[str_detect(tmp,"bed")]
#tmp1 <- paste0(seqid, ".bed")
#seqid[!(tmp1 %in% tmp)]
## [1] "SeqId_13513_174" "SeqId_16852_10"  "SeqId_5981_6"    "SeqId_7967_38"  "SeqId_6294_11"
#writeLines(seqid[(tmp1 %in% tmp)], "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt")
#
#

tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq_remove_ambiguous_snp")
library(stringr)
length(tmp[str_detect(tmp,"bed")]) #[1] 4645
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
length(seqid) #[1] 4665
tmp <- tmp[str_detect(tmp,"bed")]
tmp1 <- paste0(seqid, ".bed")
seqid[!(tmp1 %in% tmp)]
 #[1] "SeqId_11682_7"   "SeqId_13513_174" "SeqId_14204_55"  "SeqId_15482_12"
 #[5] "SeqId_16852_10"  "SeqId_5981_6"    "SeqId_6440_31"   "SeqId_6911_103"
 #[9] "SeqId_7752_31"   "SeqId_7967_38"   "SeqId_8853_2"    "SeqId_6294_11"
writeLines(seqid[(tmp1 %in% tmp)], "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt")

#
#
############################################################################
############################################################################
### cis-SNPs
#
#seqid <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt")
#
#res <- 0
#for (i in 1:length(seqid)){
#  res <- res + dim(read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq_remove_ambiguous_snp/",seqid[i],".bim")))[1]
#  print(i)
#}
#res
#
#
