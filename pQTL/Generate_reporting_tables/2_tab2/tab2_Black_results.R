
####################################
####################################
## conditional report

## Black

rm(list=ls())

library(readr)
ethnic <- "Black"

tab2  <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/tab2.txt"))

SOMA <- unique(tab2$SOMAmer)
length(SOMA)

R2 <- numeric()
for (i in 1:length(SOMA)){
  a <- readRDS( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/tab2/R2_", i, ".rds"))
  R2 <- c(R2, a)
  print(i)
}

tab2 <- cbind(tab2, R2)
write_tsv(tab2,"/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/tab2/tab2_R2.txt")