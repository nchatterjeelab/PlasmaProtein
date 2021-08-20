library(readr)
library(stringr)
b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/Plasma_Protein_AA_hg38.pos")
tmp <- paste0("chr",b$CHR,":",b$P0,"-",b$P0)
writeLines(tmp, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg38.bed")
tmp0 <- b$P0
tmp0[tmp != "chr1:206008535-206008535"] <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg37_AA_P0.bed")
tmp0[tmp == "chr1:206008535-206008535"] <- NA
b$P0 <- unlist(lapply(str_split(tmp0,":|-"), function(x){x[2]}))

tmp <- paste0("chr",b$CHR,":",b$P1,"-",b$P1)
writeLines(tmp, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg38.bed")
tmp0 <- b$P1
tmp0[tmp != "chr7:142774564-142774564"] <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg37_AA_P1.bed")
tmp0[tmp == "chr7:142774564-142774564"] <- NA
b$P1 <- unlist(lapply(str_split(tmp0,":|-"), function(x){x[2]}))

write_tsv(b,"/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/Plasma_Protein_AA_hg19.pos")


library(readr)
library(stringr)
b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/Plasma_Protein_EA_hg38.pos")
tmp <- paste0("chr",b$CHR,":",b$P0,"-",b$P0)
writeLines(tmp, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg38.bed")
tmp0 <- b$P0
tmp0[tmp != "chr1:206008535-206008535"] <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg37_EA_P0.bed")
tmp0[tmp == "chr1:206008535-206008535"] <- NA
b$P0 <- unlist(lapply(str_split(tmp0,":|-"), function(x){x[2]}))

tmp <- paste0("chr",b$CHR,":",b$P1,"-",b$P1)
writeLines(tmp, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg38.bed")
tmp0 <- b$P1
tmp0[tmp != "chr7:142774564-142774564"] <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/hg37_EA_P1.bed")
tmp0[tmp == "chr7:142774564-142774564"] <- NA
b$P1 <- unlist(lapply(str_split(tmp0,":|-"), function(x){x[2]}))

write_tsv(b,"/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/Plasma_Protein_EA_hg19.pos")
