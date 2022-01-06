
## Fig 3

library(stringr)
library(readr)
library(dplyr)
library(xlsx)


summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")


###############################################################
###############################################################
###############################################################

#### A

summary.Liver <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Liver.rds")
summary.Whole_Blood <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Whole_Blood.rds")

summary.Whole_Blood$M <- length(summary.Whole_Blood$RDat)
summary.Liver$M <- length(summary.Liver$RDat)
df.hsq <- data.frame(hsq = c(summary.AA$hsq.all, summary.EA$hsq.all,
                             summary.Liver$hsq.all, summary.Whole_Blood$hsq.all), 
                     group=c(rep("P in AA", summary.AA$M), rep("P in EA", summary.EA$M),
                             rep("T in Liver", summary.Liver$M), rep("T in Whole Blood", summary.Whole_Blood$M)),
                     kind=c(rep("Plasma protein\nSOMAmers (P)", summary.AA$M), rep("Plasma protein\nSOMAmers (P)", summary.EA$M),
                            rep("Gene expr-\nessions (T)", summary.Liver$M), rep("Gene expr-\nessions (T)", summary.Whole_Blood$M)),
                     stringsAsFactors = F)
df.hsq$group <- factor(df.hsq$group, levels = c("T in Liver", "T in Whole Blood", "P in AA", "P in EA"))

write.xlsx(df.hsq, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig3.xlsx",
           sheetName = "3a", row.names = FALSE, append = FALSE)


###############################################################
###############################################################
###############################################################

##### B

library(RColorBrewer)

df.acc <- data.frame(acc = c(summary.AA$rsq.enet/summary.AA$hsq.all, summary.EA$rsq.enet/summary.EA$hsq.all,
                             summary.Liver$rsq.enet/summary.Liver$hsq.all, summary.Whole_Blood$rsq.enet/summary.Whole_Blood$hsq.all,
                             summary.AA$rsq.top1/summary.AA$hsq.all, summary.EA$rsq.top1/summary.EA$hsq.all,
                             summary.Liver$rsq.top1/summary.Liver$hsq.all, summary.Whole_Blood$rsq.top1/summary.Whole_Blood$hsq.all), 
                     group=c(rep("P in AA", summary.AA$M), rep("P in EA", summary.EA$M),
                             rep("T in Liver", summary.Liver$M), rep("T in Whole Blood", summary.Whole_Blood$M),
                             rep("P in AA", summary.AA$M), rep("P in EA", summary.EA$M),
                             rep("T in Liver", summary.Liver$M), rep("T in Whole Blood", summary.Whole_Blood$M)),
                     Model=c(rep("Elastic Net", summary.AA$M+summary.EA$M+summary.Liver$M+summary.Whole_Blood$M),
                             rep("Best cis-SNP", summary.AA$M+summary.EA$M+summary.Liver$M+summary.Whole_Blood$M)), 
                     stringsAsFactors = F)
df.acc$group <- factor(df.acc$group, levels = c("T in Liver", "T in Whole Blood", "P in AA", "P in EA"))
df.acc$Model <- factor(df.acc$Model, levels=c("Best cis-SNP", "Elastic Net"))

write.xlsx(df.acc, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig3.xlsx",
           sheetName = "3b", row.names = FALSE, append = TRUE)


###############################################################
###############################################################
###############################################################

#### C

pred_acc <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/cross_ethnic/EAtoAA_Predicted_Accuracy.txt")
pred_acc$acc <- pred_acc$rsq/pred_acc$hsq
pred_acc$acc.cross <- pred_acc$rsq.cross/pred_acc$hsq
pred_acc$ratio <- pred_acc$rsq.cross/pred_acc$rsq
df.EAtoAA <- data.frame(acc = c(pred_acc$acc, pred_acc$acc.cross),
                        model=c(rep("AA Model", dim(pred_acc)[1]), 
                                rep("EA Model", dim(pred_acc)[1] )))

pred_acc <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/cross_ethnic/AAtoEA_Predicted_Accuracy.txt")
pred_acc$acc <- pred_acc$rsq/pred_acc$hsq
pred_acc$acc.cross <- pred_acc$rsq.cross/pred_acc$hsq
pred_acc$ratio <- pred_acc$rsq.cross/pred_acc$rsq
df.AAtoEA <- data.frame(acc = c(pred_acc$acc, pred_acc$acc.cross),
                        model=c(rep("EA Model", dim(pred_acc)[1]), 
                                rep("AA Model", dim(pred_acc)[1] )))

df.CE <- cbind(rbind(df.EAtoAA, df.AAtoEA), 
               race=c(rep("AA",nrow(df.EAtoAA)),rep("EA",nrow(df.AAtoEA))))

write.xlsx(df.CE, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig3.xlsx",
           sheetName = "3c", row.names = FALSE, append = TRUE)


###############################################################
###############################################################
###############################################################

#### D

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)

res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/correlations_v7.txt")

a <- unlist(lapply(strsplit(res$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
res$tissue <- gtex.colors$V1[match(a, b)]
tmp <- res %>% group_by(tissue) %>% summarise(med = median(cor, na.rm=T))
tmp <- tmp[order(tmp$med,decreasing=T),]
res$tissue = factor(res$tissue, levels=tmp$tissue)
res <- as.data.frame(res)
write.xlsx(res, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig3.xlsx",
           sheetName = "3d", row.names = FALSE, append = TRUE)


