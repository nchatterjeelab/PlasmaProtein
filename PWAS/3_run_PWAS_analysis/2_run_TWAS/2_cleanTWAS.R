
########################################################################
## V7

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7.out"))
tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V7_tissue_list.txt")

library(dplyr)
library(readr)

for (tissue in tissue_list){

  results <- tibble()
  for (chr in 1:22) {
      results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7/",tissue,"/chr", chr, ".out")))
      if(chr==6){
          results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7/",tissue,"/chr", chr, ".out.MHC")))
      }
  }
  colnames(results) <- c("PANEL","FILE","ID","CHR","P0","P1","HSQ",
                         "BEST.GWAS.ID","BEST.GWAS.Z",
                         "EQTL.ID","EQTL.R2","EQTL.Z","EQTL.GWAS.Z",
                         "NSNP","NWGT","MODEL","MODELCV.R2","MODELCV.PV",
                         "TWAS.Z","TWAS.P","TWAS.BETA","TWAS.SE","TWAS.CI")

  write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7.out/",tissue,".out"))

}


#### merge all tissues (TWAS)
#
#tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_tissue_list.txt")
#library(dplyr)
#library(readr)
#library(stringr)
#results <- tibble()
#for (tissue in tissue_list){
#
#  results <- rbind(results,read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7.out/",tissue,".out"), col_types = "ccciiidcdcdddiicddddddc"))
#  print(tissue)
#
#}
#write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7.out/all-cleaned.out"))


########################################################################
## V8

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8.out"))
tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_tissue_list.txt")

library(dplyr)
library(readr)

for (tissue in tissue_list){

  results <- tibble()
  for (chr in 1:22) {
      results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8/",tissue,"/chr", chr, ".out")))
      if(chr==6){
          results <- rbind(results, read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8/",tissue,"/chr", chr, ".out.MHC")))
      }
  }
  colnames(results) <- c("PANEL","FILE","ID","CHR","P0","P1","HSQ",
                         "BEST.GWAS.ID","BEST.GWAS.Z",
                         "EQTL.ID","EQTL.R2","EQTL.Z","EQTL.GWAS.Z",
                         "NSNP","NWGT","MODEL","MODELCV.R2","MODELCV.PV",
                         "TWAS.Z","TWAS.P","TWAS.BETA","TWAS.SE","TWAS.CI")

  write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8.out/",tissue,".out"))

}


#### merge all tissues (TWAS)
#
#tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_tissue_list.txt")
#library(dplyr)
#library(readr)
#library(stringr)
#results <- tibble()
#for (tissue in tissue_list){
#
#  results <- rbind(results,read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8.out/",tissue,".out"), col_types = "ccciiidcdcdddiicddddddc"))
#  print(tissue)
#
#}
#write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8.out/all-cleaned.out"))
#
