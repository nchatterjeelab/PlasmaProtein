
## Fig 4

library(readr)
library(xlsx)

###############################################################
###############################################################
###############################################################

# Urate (panel a)

disease <- "Urate"

#####################
## data preparation

dat.pwas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_",disease,"_CI_hg19.txt"))
dat.pwas <- dat.pwas[!(is.na(dat.pwas$PWAS.P)),]
dat.pwas <- dat.pwas[!(is.na(dat.pwas$P0)),]
dat.pwas <- dat.pwas
dat.pwas$PWAS.P[dat.pwas$ID=="INHBC"] = 10^(-30)
dat.pwas$ID[dat.pwas$ID=="INHBC"] = "INHBC(7.95e-63)"

dat.pwas <- dat.pwas[,c("ID","CHR","P0","PWAS.P")]
colnames(dat.pwas) <- c("ID","CHR","BP","P")
dat.pwas$tissue <- rep("Plasma",nrow(dat.pwas))

dat.pwas <- dat.pwas[!(is.na(dat.pwas$P)),]
dat.pwas$pthres <- 0.05/nrow(dat.pwas)


library(dplyr)
tissue_list <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/Pipeline/PWAS/GTEx_V7_tissue_list.txt")
res <- tibble()
for(tissue in tissue_list){
  dat.twas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/TWAS/", disease, "_EA/", tissue, ".out"))
  dat.twas <- dat.twas[!(is.na(dat.twas$TWAS.P)),]
  dat.twas$tissue <- rep(tissue,nrow(dat.twas))
  res <- rbind(res,dat.twas)
}

dat.twas <- res[,c("ID","CHR","P0","TWAS.P","tissue")]
colnames(dat.twas) <- c("ID","CHR","BP","P","tissue")

dat.twas <- dat.twas[!(is.na(dat.twas$P)),]
dat.twas$pthres <- 0.05/nrow(dat.twas)

dat_all <- rbind(dat.pwas, dat.twas)


nCHR <- length(unique(dat_all$CHR))
dat_all$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(dat_all$CHR)){
  nbp[i] <- max(dat_all[dat_all$CHR == i,]$BP)
  dat_all$BPcum[dat_all$CHR == i] <- dat_all$BP[dat_all$CHR == i] + s
  s <- s + nbp[i]
}
dat_all$disease <- disease

dat_all_disease <- dat_all

###############################################################
###############################################################
###############################################################

# Gout (panel b)


disease <- "Gout"

#####################
## data preparation

dat.pwas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_",disease,"_CI_hg19.txt"))
dat.pwas <- dat.pwas[!(is.na(dat.pwas$PWAS.P)),]
dat.pwas <- dat.pwas[!(is.na(dat.pwas$P0)),]
dat.pwas <- dat.pwas

dat.pwas <- dat.pwas[,c("ID","CHR","P0","PWAS.P")]
colnames(dat.pwas) <- c("ID","CHR","BP","P")
dat.pwas$tissue <- rep("Plasma",nrow(dat.pwas))

dat.pwas <- dat.pwas[!(is.na(dat.pwas$P)),]
dat.pwas$pthres <- 0.05/nrow(dat.pwas)


library(dplyr)
tissue_list <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/Pipeline/PWAS/GTEx_V7_tissue_list.txt")
res <- tibble()
for(tissue in tissue_list){
  dat.twas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/TWAS/", disease, "_EA/", tissue, ".out"))
  dat.twas <- dat.twas[!(is.na(dat.twas$TWAS.P)),]
  dat.twas$tissue <- rep(tissue,nrow(dat.twas))
  res <- rbind(res,dat.twas)
}

dat.twas <- res[,c("ID","CHR","P0","TWAS.P","tissue")]
colnames(dat.twas) <- c("ID","CHR","BP","P","tissue")

dat.twas <- dat.twas[!(is.na(dat.twas$P)),]
dat.twas$pthres <- 0.05/nrow(dat.twas)

dat_all <- rbind(dat.pwas, dat.twas)



nCHR <- length(unique(dat_all$CHR))
dat_all$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(dat_all$CHR)){
  nbp[i] <- max(dat_all[dat_all$CHR == i,]$BP)
  dat_all$BPcum[dat_all$CHR == i] <- dat_all$BP[dat_all$CHR == i] + s
  s <- s + nbp[i]
}
dat_all$disease <- disease

dat_all_disease <- rbind(dat_all_disease, dat_all)

write_tsv(dat_all_disease, "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig4.txt")


