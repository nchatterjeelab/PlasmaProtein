
rm(list=ls())

#library(readr)
#library(plink2R)
library(bigreadr)
#library(stringr)


aask <- character()
for (chr in 1:22){
  aask <- c(aask, fread2(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink/chr",chr,".bim"))$V2)
  print(chr)
}
writeLines(aask,paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink/all_snp_ID.txt"))

aric_b <- character()
for (chr in 1:22){
  aric_b <- c(aric_b, fread2(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/chr",chr,".bim"))$V2)
  print(chr)
}
writeLines(aric_b,paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/Black/all_snp_ID.txt"))


aric_w <- character()
for (chr in 1:22){
  aric_w <- c(aric_w, fread2(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/chr",chr,".bim"))$V2)
  print(chr)
}
writeLines(aric_w,paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/all_snp_ID.txt"))

SNP <- unique(c(aask,aric_b,aric_w))
writeLines(SNP, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/plink/all_snp_ID_aaskandaric.txt"))

lookup <- fread2("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/rsid/GRCh38_dbSNP151_rsid_final_USE.txt")
a <- lookup[lookup$"chr1:10531:C:G" %in% SNP,]
colnames(a) <- c("SNPid","rsid")
saveRDS(a, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")
write_tsv(a, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.txt", col_names=F)

a <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")

b <- table(a$SNPid)
c <- names(b)[b>1]
readr::write_tsv(a[a$SNPid %in% c,], "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_dup.txt", col_names=F)
writeLines(c,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_dup_id.txt")

a <- a[!(a$SNPid %in% c),]
saveRDS(a, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_nodup.rds")
readr::write_tsv(a, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_nodup.txt", col_names=F)




