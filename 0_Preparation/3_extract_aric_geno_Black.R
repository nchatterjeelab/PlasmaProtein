
################################################################
################################################################
################################################################
# 1M window

rm(list=ls())

library(readr)

for (chr in 1:22){
 a <- read.table(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr",chr,".bim"), stringsAsFactors = F)
 ref1 <- a$V5
 flip1 = ref1
 flip1[ref1 == "A"] = "T"
 flip1[ref1 == "T"] = "A"
 flip1[ref1 == "G"] = "C"
 flip1[ref1 == "C"] = "G"

 ref2 <- a$V6
 writeLines(a$V2[!(ref2==flip1)],paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Remove_strand_ambiguous/Black/Remaining_SNP/chr",chr))
 print(chr)
}

for (chr in 1:22){
 system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 --bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr", chr, " --extract /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Remove_strand_ambiguous/Black/Remaining_SNP/chr",chr," --make-bed --out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Remove_strand_ambiguous/Black/chr",chr))
}

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')


chr <- annota$chromosome_name
TSS <- annota$transcription_start_site

dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M")
dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp")

b <- "#!/usr/bin/env bash
#$ -N window1M
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0(b,
"

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Remove_strand_ambiguous/Black/chr",chr[i], " \\
--extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_GRCh38_ARICID.txt \\
--chr ", chr[i], " --from-bp ", TSS[i]-500000," --to-bp ", TSS[i]+500000, " \\
--force-intersect \\
--make-bed \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", seqid[i], "

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", seqid[i], " \\
--update-name /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_to_RSID.txt \\
--make-bed \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", seqid[i], "

rm /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", seqid[i], ".bed~
rm /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", seqid[i], ".bim~
rm /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", seqid[i], ".fam~

####################

")
 print(i)
}

writeLines(b,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/extractgene_remove_ambiguous_snp.sh')


tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp")
library(stringr)
length(tmp[str_detect(tmp,"bed")]) #[1] 4645
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
length(seqid) #[1] 4665
tmp <- tmp[str_detect(tmp,"bed")]
tmp1 <- paste0(seqid, ".bed")
seqid[!(tmp1 %in% tmp)]
# [1] "SeqId_11682_7"   "SeqId_13513_174" "SeqId_14204_55"  "SeqId_15482_12"
# [5] "SeqId_16755_195" "SeqId_16852_10"  "SeqId_2732_58"   "SeqId_3332_57"
# [9] "SeqId_3708_62"   "SeqId_5028_59"   "SeqId_5981_6"    "SeqId_6440_31"
#[13] "SeqId_6580_29"   "SeqId_6911_103"  "SeqId_7752_31"   "SeqId_7960_53"
#[17] "SeqId_7967_38"   "SeqId_8853_2"    "SeqId_9094_5"    "SeqId_6294_11"
writeLines(seqid[(tmp1 %in% tmp)], "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/seqid_autosomal_withSNP.txt")



###########################################################################
###########################################################################

# 2M window


rm(list=ls())

library(readr)


annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')


chr <- annota$chromosome_name
TSS <- annota$transcription_start_site

dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M")
dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq_remove_ambiguous_snp")

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
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr", chr[i], " \\
--chr ", chr[i], " --from-bp ", TSS[i]-1000000," --to-bp ", TSS[i]+1000000, " \\
--make-bed \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq_remove_ambiguous_snp/", seqid[i], "

####################

")
print(i)
}

writeLines(b,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/extractgene.sh')




tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/byseq_remove_ambiguous_snp")
library(stringr)
length(tmp[str_detect(tmp,"bed")]) #[1] 4660
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
length(seqid) #[1] 4665
tmp <- tmp[str_detect(tmp,"bed")]
tmp1 <- paste0(seqid, ".bed")
seqid[!(tmp1 %in% tmp)]
# [1] "SeqId_13513_174" "SeqId_16852_10"  "SeqId_5981_6"    "SeqId_7967_38"  "SeqId_6294_11"
writeLines(seqid[(tmp1 %in% tmp)], "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window2M/seqid_autosomal_withSNP.txt")



