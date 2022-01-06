
## Suppl Fig 8

library(readr)

dat <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_Urate_CI.txt"))
dat1 <- dat[!is.na(dat$PWAS.Z),]
dat1$disease <- "Urate"

dat <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_Gout_CI.txt"))
dat2 <- dat[!is.na(dat$PWAS.Z),]
dat2$disease <- "Gout"

dat_disease <- rbind(dat1, dat2)
write_tsv(dat_disease, paste0("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig8.txt"))


