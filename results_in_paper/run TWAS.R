########################################################################

# TWAS


## Urate_EA summary data for TWAS
library(bigreadr)

dat <- read.table("/dcl01/chatterj/data/jzhang2/summary_data/Urate/EA.txt", header = T, sep=" ")
dat1 <- dat[,c(3,4,5,7,8)]
dat1$Z <- dat1$Effect/dat1$StdErr
colnames(dat1) <- c("SNP", "A1", "A2","BETA","SE","Z")
tmp <- as.character(dat1$SNP)
dat1 <- dat1[tmp!="", ]
readr::write_tsv(dat1, "/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Urate_EA_cleaned_rsid_BETA_SE_TWAS.txt")

system("cp /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Urate_EA_cleaned_rsid_BETA_SE_TWAS.txt /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Urate_EA_cleaned_rsid_BETA_SE_PWAS.txt")


## Gout_EA summary data for TWAS

rm(list=ls())

library(readr)
library(dplyr)
library(stringr)

dat <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/summary_data/Gout/gout_chr1_22_LQ_IQ06_mac10_EA_171.tbl.nstud9.summac400"))
snp.ref <- read_tsv("/dcl01/chatterj/data/jzhang2/TWAS/snp151Common_cleaned.txt")

summ=dat

a <- str_split(summ$MarkerName,"_")

summ$Chr <- as.integer(unlist(lapply(a,FUN=function (x){x[1]})))
summ$Pos <- as.integer(unlist(lapply(a,FUN=function (x){x[2]})))

summ <- summ[,c("Chr", "Pos", "Allele1", "Allele2", "Effect", "StdErr")]
summ$snp.key <- paste0("chr",summ$Chr,":",format(summ$Pos, trim = T, scientific = F))

summ <- inner_join(summ, snp.ref[,4:5],by="snp.key")
summ$Z <- summ$Effect/summ$StdErr
summ <- summ[,c("rsid", "Allele1", "Allele2", "Z","Effect", "StdErr")]

colnames(summ) <- c("SNP", "A1", "A2", "Z","BETA","SE")

write_tsv(summ, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Gout_EA_cleaned_rsid_BETA_SE_TWAS.txt"))


######################################################################
## run TWAS

rm(list=ls())
dir.create("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/submit/Gout/TWAS_CI")

tissue_list <- readLines("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/GTEx_V7_tissue_list.txt")

b <- paste0("#!/usr/bin/env bash
#$ -N Gout
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (tissue in tissue_list){
dir.create(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI/",tissue))

a <- paste0("#!/usr/bin/env bash
#$ -N Gout_",tissue,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript /dcl01/chatterj/data/jzhang2/PWAS_tutorial/scripts/PWAS.assoc_test_CI.R \\
--sumstats /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Gout_EA_cleaned_rsid_BETA_SE_TWAS.txt \\
--weights /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/",tissue,".P01.pos \\
--weights_dir /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/ \\
--ref_ld_chr /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/LDREF/1000G.EUR. \\
--force_model enet \\
--chr $SGE_TASK_ID \\
--out /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI/",tissue,"/chr$SGE_TASK_ID.out

")
writeLines(a, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/submit/Gout/TWAS_CI/",tissue,".sh"))

  b <- paste0(b, "
qsub ",tissue,".sh")

}
writeLines(b, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/submit/Gout/TWAS_CI/ALL.sh"))


########################################################################
## clean table

dir.create(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI.out"))
tissue_list <- readLines("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/GTEx_V7_tissue_list.txt")
library(dplyr)
library(readr)

for (tissue in tissue_list){

  results <- tibble()
  for (chr in 1:22) {
      results <- rbind(results, read_tsv(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI/",tissue,"/chr", chr, ".out")))
      if(chr==6){
          results <- rbind(results, read_tsv(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI/",tissue,"/chr", chr, ".out.MHC")))
      }
  }
  colnames(results) <- c("PANEL","FILE","ID","CHR","P0","P1","HSQ",
                         "BEST.GWAS.ID","BEST.GWAS.Z",
                         "EQTL.ID","EQTL.R2","EQTL.Z","EQTL.GWAS.Z",
                         "NSNP","NWGT","MODEL","MODELCV.R2","MODELCV.PV",
                         "TWAS.Z","TWAS.P","TWAS.BETA","TWAS.SE","TWAS.CI")

  write_tsv(results, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI.out/",tissue,".out"))

}


### merge all tissues


tissue_list <- readLines("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/GTEx_V7_tissue_list.txt")
library(dplyr)
library(readr)
library(stringr)
results <- tibble()
for (tissue in tissue_list){

  results <- rbind(results,read_tsv(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI.out/",tissue,".out"), col_types = "ccciiidcdcdddiicddddddc"))
  print(tissue)

}
write_tsv(results, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI.out/all-cleaned.out"))



########################################################################
## conditional analysis

"
Rscript /dcl01/chatterj/data/jzhang2/PWAS_tutorial/scripts/PWAS.conditional_CI.R \
--PWAS /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/PWAS_CI.out \
--TWAS /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/TWAS_CI.out/ \
--tissue_list /dcl01/chatterj/data/jzhang2/PWAS_tutorial/GTEx_V7_tissue_list.txt \
--tissue_n_gene /dcl01/chatterj/data/jzhang2/PWAS_tutorial/GTEx_V7_n_gene.rds \
--imputed_P /dcl01/chatterj/data/jzhang2/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_Plasma_Protein.txt \
--imputed_T /dcl01/chatterj/data/jzhang2/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_FUSION/ \
--out /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Gout/ConditionalAnalysis_CI/
"


