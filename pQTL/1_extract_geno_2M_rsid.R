
#for (chr in 1:22){
# system(paste0("/users/jzhang2/RESEARCH/tools/plink/plink2",
#               " --threads 1",
#               " --bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Remove_strand_ambiguous/White/chr", chr,
#               " --update-name /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_nodup.txt",
#               " --make-bed",
#               " --out /fastscratch/myscratch/jzhang2/ARIC/chr", chr))
#}
#

rm(list=ls())

library(readr)

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')

annota <- annota[match(seqid,annota$seqid_in_sample),]

chr <- annota$chromosome_name
TSS <- annota$transcription_start_site

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/byseq_remove_ambiguous_snp")


b <- "#!/usr/bin/env bash
#$ -N window2M
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
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Remove_strand_ambiguous/White/update_rsid/chr", chr[i], " \\
--chr ", chr[i], " --from-bp ", TSS[i]-1000000," --to-bp ", TSS[i]+1000000, " \\
--make-bed \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/byseq_remove_ambiguous_snp/", seqid[i], "

")
 print(i)
}

writeLines(b,  '/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/extractgene_remove_ambiguous_snp.sh')





tmp <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/byseq_remove_ambiguous_snp")
library(stringr)
length(tmp[str_detect(tmp,"bed")]) #[1] 4657
tmp <- tmp[str_detect(tmp,"bed")]
tmp <- gsub(".bed","",tmp)
writeLines(tmp, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M/seqid_autosomal_withSNP_remove_ambiguous_snp.txt")

