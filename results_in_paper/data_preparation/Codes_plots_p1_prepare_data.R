
## Fig 1

#### A

library(readxl)
library(xlsx)

tmp <- read_excel("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Suppl_tables.xlsx", sheet = "ST2 -- PEER factors", skip=3)
tmp <- tmp[1:20,]
df.peer <- data.frame(NPEER=rep(tmp[[1]],2), NpGene=c(tmp[[2]], tmp[[3]]), Race=c(rep("EA", nrow(tmp)), rep("AA", nrow(tmp))))
df.peer$NPEER <- as.integer(as.character(df.peer$NPEER))
df.peer <- df.peer[df.peer$NPEER<140,]

write.xlsx(df.peer, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig1.xlsx",
           sheetName = "1a", row.names = FALSE, append = FALSE)


###############################################################
###############################################################
###############################################################

##### B

df.venn <- data.frame(x = c(-0.1,0.1),
                      y = c(0,0),
                      r = c(0.4,0.45),
                      labels = c('A', 'B'))
write.xlsx(df.venn, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig1.xlsx",
           sheetName = "1b", row.names = FALSE, append = TRUE)

###############################################################
###############################################################
###############################################################

## C

sig.w <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_pQTL_summary_cleaned_4.0_EA.txt")
sig.b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_pQTL_summary_cleaned_4.0_AA.txt")
library(dplyr)
tmp <- read.table("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_permutations_all.thresholds_EA.txt")
sig.w <- inner_join(sig.w, tmp, by=c("SOMAmer"="V1"))
sig.w$freqvar <- sig.w$A1_AF*(1-sig.w$A1_AF)
sig.w$fR2 <- 2 * sig.w$freqvar * sig.w$Beta^2
sig.w$power <- pf(qf(sig.w$V2, df1 = 1, df2 = 7213-1-1, lower.tail = FALSE),
                  df1 = 1, df2 = 7213-1-1, lower.tail = FALSE,
                  ncp = 7213*sig.w$fR2/(1-sig.w$fR2))
tmp <- read.table("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_permutations_all.thresholds_AA.txt")
sig.b <- inner_join(sig.b, tmp, by=c("SOMAmer"="V1"))
sig.b$freqvar <- sig.b$A1_AF*(1-sig.b$A1_AF)
sig.b$fR2 <- 2 * sig.b$freqvar * sig.b$Beta^2
sig.b$power <- pf(qf(sig.b$V2, df1 = 1, df2 = 1871-1-1, lower.tail = FALSE),
                  df1 = 1, df2 = 1871-1-1, lower.tail = FALSE,
                  ncp = 1871*sig.b$fR2/(1-sig.b$fR2))

df.effmaf <- data.frame(BETA = c(sig.w$Beta, sig.b$Beta), 
                        maf = c(sig.w$freqvar, sig.b$freqvar),
                        power = c(sig.w$power, sig.b$power),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))

write.xlsx(df.effmaf, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig1.xlsx",
           sheetName = "1c", row.names = FALSE, append = TRUE)


###############################################################
###############################################################
###############################################################

## D

df.efftss <- data.frame(BETA = c(sig.w$Beta, sig.b$Beta), 
                        dist = c(sig.w$TopSNP_Pos38-sig.w$TSS, sig.b$TopSNP_Pos38-sig.b$TSS),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))

write.xlsx(df.efftss, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig1.xlsx",
           sheetName = "1d", row.names = FALSE, append = TRUE)


###############################################################
###############################################################
###############################################################

## E

sig.w <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/2_conditional_analysis_1.0_EA.txt")
sig.b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/2_conditional_analysis_1.0_AA.txt")


df.cond <- data.frame(count=c(data.frame(sort(table(sig.w$SOMAmer),decreasing=T))$Freq, 
                              data.frame(sort(table(sig.b$SOMAmer),decreasing=T))$Freq),
                      Race=c(rep("EA", length(data.frame(sort(table(sig.w$SOMAmer),decreasing=T))$Freq)),
                             rep("AA", length(data.frame(sort(table(sig.b$SOMAmer),decreasing=T))$Freq))))

write.xlsx(df.cond, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig1.xlsx",
           sheetName = "1e", row.names = FALSE, append = TRUE)


