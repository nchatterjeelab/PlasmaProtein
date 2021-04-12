
a <- dir("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_EA")
a <- a[stringr::str_detect(a,"RDat")]
all <- character()
for (i in 1:length(a)){
  load(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_EA/",a[i]))
  all <- c(all,snps$V2)
  all <- unique(all)
  print(i)
}
writeLines(all, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/all_snps.txt")


########################################

for (chr in 1:22){

  tmp <- paste0("/users/jzhang2/RESEARCH/tools/plink/plink2",
                " --extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/all_snps.txt",
                " --keep /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/ids.list",
                " --bfile /dcl01/chatterj/data/ukbiobank_extracted/imputed_bed/ukb_imp_chr",chr,"_v3",
                " --snps-only",
                " --rm-dup exclude-all",
                " --pheno /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/cont_normal.pheno",
                " --covar /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/covs.list",
                " --glm cols=+a1freq,-a1countcc hide-covar --vif 10000 --covar-variance-standardize",
                " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/cont_normal_",chr)

tmp <- paste0("#!/usr/bin/env bash
#$ -N gwas_", chr,"
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -pe local 5
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

", tmp)

  writeLines(tmp, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS_on_null_phenotype/cont_normal/chr",chr,".sh"))

}



for (chr in 1:22){

  tmp <- paste0("/users/jzhang2/RESEARCH/tools/plink/plink2",
                " --extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/all_snps.txt",
                " --keep /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/ids.list",
                " --bfile /dcl01/chatterj/data/ukbiobank_extracted/imputed_bed/ukb_imp_chr",chr,"_v3",
                " --snps-only",
                " --rm-dup exclude-all",
                " --pheno /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/binary_equal.pheno",
                " --covar /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/covs.list",
                " --glm cols=+a1freq,-a1countcc hide-covar --vif 10000 --covar-variance-standardize",
                " --1",
                " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/binary_equal_",chr)

tmp <- paste0("#!/usr/bin/env bash
#$ -N gwas_", chr,"
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 5
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

", tmp)

  writeLines(tmp, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS_on_null_phenotype/binary_equal/chr",chr,".sh"))

}



for (chr in 1:22){

  tmp <- paste0("/users/jzhang2/RESEARCH/tools/plink/plink2",
                " --extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/all_snps.txt",
                " --keep /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/ids.list",
                " --bfile /dcl01/chatterj/data/ukbiobank_extracted/imputed_bed/ukb_imp_chr",chr,"_v3",
                " --snps-only",
                " --rm-dup exclude-all",
                " --pheno /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/binary_10x.pheno",
                " --covar /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/covs.list",
                " --glm cols=+a1freq,-a1countcc hide-covar --vif 10000 --covar-variance-standardize",
                " --1",
                " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/binary_10x_",chr)

tmp <- paste0("#!/usr/bin/env bash
#$ -N gwas_", chr,"
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 5
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

", tmp)

  writeLines(tmp, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS_on_null_phenotype/binary_10x/chr",chr,".sh"))

}

library(readr)
for (chr in 1:22){
  a <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/cont_normal_",chr,".PHENO1.glm.linear"))
  a <- a[,c(3:6,10,11)]
  a$A2 <- a$ALT
  tmp <- a$REF != a$A1
  a$A2[tmp] <- a$REF[tmp]
  a$Z <- a$BETA/a$SE
  a <- a[,c(1,4,7,8)]
  colnames(a) <- c("SNP","A1","A2","Z")
  write_tsv(a,paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/cont_normal_",chr,"_cleaned"))
}

library(dplyr)
results <- tibble()
for (chr in 1:22) {
    results <- rbind(results, read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/PWAS/cont_normal_chr",chr,".out")))
    if(chr==6){
        results <- rbind(results, read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/PWAS/cont_normal_chr",chr,".out.MHC")))
    }
}
write_tsv(results, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/PWAS/cont_normal_allchr.out"))

