

##################################################
## Clean existing studies data

library(readxl)
library(dplyr)
library(stringr)
library(readr)

####################

pQTL <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_pQTL_summary_cleaned_4.0_EA.txt")
coloc <- read_excel("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_06_revision2/Suppl_tables_9Aug2021_DD_JZ.xlsx",
                   sheet = "ST8.2-- coloc(PP.H4) with eQTLs", skip = 2)
m <- as.matrix(coloc[,-1:-2])
class(m) <- "numeric"

pph4 <- apply(m, MARGIN = 1, function(x){max(x, na.rm=T)})
pph4[pph4==-Inf] <- NA
pQTL$pph4 <- pph4
write_tsv(pQTL, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_pQTL_summary_cleaned_5.0_EA.txt")


