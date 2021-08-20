
a <- dir("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA")
a <- a[stringr::str_detect(a,"RDat")]
all <- character()
for (i in 1:length(a)){
  load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA/",a[i]))
  all <- c(all,snps$V2)
  all <- unique(all)
  print(i)
}
writeLines(all, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/all_snps.txt")


########################################

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/cont_normal/")
tmp <- paste0("#!/usr/bin/env bash
#$ -N gwas_continuous
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -pe local 10
#$ -t 1-22
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 10",
                " --extract /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/all_snps.txt",
                " --keep /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/ids.list",
                " --bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3",
                " --snps-only",
                " --rm-dup exclude-all",
                " --pheno /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/cont_normal.pheno",
                " --covar /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/covs.list",
                " --glm cols=+a1freq,-a1countcc hide-covar --vif 10000 --covar-variance-standardize",
                " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/cont_normal_${SGE_TASK_ID}")

writeLines(tmp, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/cont_normal/gwas.sh"))


dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/binary_equal/")
tmp <- paste0("#!/usr/bin/env bash
#$ -N gwas_binary_equal
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -pe local 10
#$ -t 1-22
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 10",
                " --extract /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/all_snps.txt",
                " --keep /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/ids.list",
                " --bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3",
                " --snps-only",
                " --rm-dup exclude-all",
                " --pheno /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/binary_equal.pheno",
                " --covar /dcl01/chatterj/data/diptavo/GEE/UKB_unrel/covs.list",
                " --glm cols=+a1freq,-a1countcc hide-covar --vif 10000 --covar-variance-standardize",
                " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/binary_equal_${SGE_TASK_ID}")

writeLines(tmp, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/binary_equal/gwas.sh"))



########################################################################
########################################################################

## clean GWAS results

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
  a <- a[!is.na(a$Z),]
  write_tsv(a,paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/GWAS/cont_normal_",chr,"_cleaned"))
}


library(readr)
for (chr in 1:22){
  a <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/nullPWAS/binary_equal_",chr,".PHENO1.glm.logistic.hybrid"))
  a <- a[,c(3:6,13)]
  a$A2 <- a$ALT
  tmp <- a$REF != a$A1
  a$A2[tmp] <- a$REF[tmp]
  a <- a[,c(1,4,6,5)]
  colnames(a) <- c("SNP","A1","A2","Z")
  a <- a[!is.na(a$Z),]
  write_tsv(a,paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/GWAS/binary_equal_",chr,"_cleaned"))
}

########################################################################
########################################################################

## run PWAS

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/cont_normal")

a <- paste0("#!/usr/bin/env bash
#$ -N PWAS_cont_normal
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/scripts/PWAS.assoc_test.R \\
--sumstats /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/GWAS/cont_normal_${SGE_TASK_ID}_cleaned \\
--weights /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos \\
--weights_dir /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA/ \\
--ref_ld_chr /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/LDref/EUR/chr \\
--force_model enet \\
--chr $SGE_TASK_ID \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/cont_normal/chr$SGE_TASK_ID.out
")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/cont_normal/PWAS.sh"))


dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/binary_equal")
a <- paste0("#!/usr/bin/env bash
#$ -N PWAS_binary_equal
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/scripts/PWAS.assoc_test.R \\
--sumstats /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/GWAS/binary_equal_${SGE_TASK_ID}_cleaned \\
--weights /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos \\
--weights_dir /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA/ \\
--ref_ld_chr /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/LDref/EUR/chr \\
--force_model enet \\
--chr $SGE_TASK_ID \\
--out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/binary_equal/chr$SGE_TASK_ID.out
")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/9_PWAS/nullPWAS/binary_equal/PWAS.sh"))


########################################################################
########################################################################

## clean PWAS results

library(dplyr)
results <- tibble()
for (chr in 1:22) {
    results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/cont_normal/chr",chr,".out")))
    if(chr==6){
        results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/cont_normal/chr",chr,".out.MHC")))
    }
}
write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/cont_normal/allchr.out"))



library(dplyr)
results <- tibble()
for (chr in 1:22) {
    results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/binary_equal/chr",chr,".out")))
    if(chr==6){
        results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/binary_equal/chr",chr,".out.MHC")))
    }
}
write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/nullPWAS/PWAS/binary_equal/allchr.out"))
