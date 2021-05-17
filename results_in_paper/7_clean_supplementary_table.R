library(readxl)
library(dplyr)

lookup <-  readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_nodup.rds")

ind <- lookup$rsid; names(ind) <- lookup$SNPid


res <- list()

dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST3.1--significant SOMAmer (EA)", skip=2)
tmp <- ind[dat$...9]; names(tmp) <- NULL
res$"ST3.1" <- tmp

dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST3.2--significant SOMAmer (AA)", skip=2)
tmp <- ind[dat$...9]; names(tmp) <- NULL
res$"ST3.2" <- tmp

dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST5.2 -- validation in AASK", skip=1)
tmp <- ind[dat$"Sentinel SNP"]; names(tmp) <- NULL
res$"ST5.2" <- tmp

dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST6.1 -- cond-indep pQTL (EA)")
tmp <- ind[dat$"sentinel SNP"]; names(tmp) <- NULL
res$"ST6.1" <- tmp
tmp <- ind[dat$"Conditional independent cis-pQTL"]; names(tmp) <- NULL
res$"ST6.1" <- cbind(res$"ST6.1",tmp)

dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST6.2 -- cond-indep pQTL (AA)")
tmp <- ind[dat$"sentinel SNP"]; names(tmp) <- NULL
res$"ST6.2" <- tmp
tmp <- ind[dat$"Conditional independent cis-pQTL"]; names(tmp) <- NULL
res$"ST6.2" <- cbind(res$"ST6.2",tmp)


dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST10.1-- single-ethnic f.m.(EA)", skip=3)
tmp <- ind[dat$"Fine-mapped SNP ID"]; names(tmp) <- NULL
res$"ST10.1" <- tmp

dat <- read_excel("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/Suppl_tables.xlsx",
                  sheet = "ST10.2-- single-ethnic f.m.(AA)", skip=3)
tmp <- ind[dat$"Fine-mapped SNP ID"]; names(tmp) <- NULL
res$"ST10.2" <- tmp

saveRDS(res, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/paper_table/rsid.rds")


################################################
################################################


library("xlsx")
res <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/*RData/rsid.rds")
tabname <- names(res)

for (i  in 1:length(tabname)) {
  write.xlsx(data.frame(res[[i]]), "/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Tables/rsid.xlsx", 
             sheetName = tabname[i], 
             col.names = T, row.names = F, append = T)
  print(i)
}
