
## Suppl Fig 5

library(readr)
library(xlsx)

###############################################################
###############################################################
###############################################################

# V7

annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/prot.anno_autosomal.txt")

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")

tmp <- summary.AA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.AA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]
tmp <- summary.EA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.EA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]

annota.Liver <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/Liver.P01.pos")
annota.Whole_Blood <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/Whole_Blood.P01.pos")

summary.Liver <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Liver.rds")
summary.Whole_Blood <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Whole_Blood.rds")

summary.Liver$gene <- annota.Liver$ID[match(paste0("Liver/",summary.Liver$RDat), annota.Liver$WGT)]
summary.Whole_Blood$gene <- annota.Whole_Blood$ID[match(paste0("Whole_Blood/",summary.Whole_Blood$RDat), annota.Whole_Blood$WGT)]

geneEALiver <- unique(intersect(summary.EA$gene, summary.Liver$gene))
geneEAWhole_Blood <- unique(intersect(summary.EA$gene, summary.Whole_Blood$gene))

geneAALiver <- unique(intersect(summary.AA$gene, summary.Liver$gene))
geneAAWhole_Blood <- unique(intersect(summary.AA$gene, summary.Whole_Blood$gene))

## AA Liver
df.hsqAALiver <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneAALiver]), 
                            group=c(rep("P in AA", sum(summary.AA$gene %in% geneAALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneAALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneAALiver)) ),
                            stringsAsFactors = F)

df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))

write.xlsx(df.hsqAALiver, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5a-1", row.names = FALSE, append = FALSE)


## AA Whole_Blood
df.hsqAAWhole_Blood <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneAAWhole_Blood]), 
                                  group=c(rep("P in AA", sum(summary.AA$gene %in% geneAAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
write.xlsx(df.hsqAAWhole_Blood, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5a-2", row.names = FALSE, append = TRUE)



## EA Liver
df.hsqEALiver <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneEALiver]), 
                            group=c(rep("P in EA", sum(summary.EA$gene %in% geneEALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneEALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneEALiver)) ),
                            stringsAsFactors = F)

df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
write.xlsx(df.hsqEALiver, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5a-3", row.names = FALSE, append = TRUE)


## EA Whole_Blood
df.hsqEAWhole_Blood <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneEAWhole_Blood]), 
                                  group=c(rep("P in EA", sum(summary.EA$gene %in% geneEAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
write.xlsx(df.hsqEAWhole_Blood, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5a-4", row.names = FALSE, append = TRUE)



###############################################################
###############################################################
###############################################################

# V8

annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/prot.anno_autosomal.txt")

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")

tmp <- summary.AA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.AA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]
tmp <- summary.EA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.EA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]

annota.Liver <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/GTEXv8.Liver.geneid.pos")
annota.Whole_Blood <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/GTEXv8.Whole_Blood.geneid.pos")

summary.Liver <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v8.ALL/summary.Liver.rds")
summary.Whole_Blood <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v8.ALL/summary.Whole_Blood.rds")

summary.Liver$gene <- annota.Liver$ID[match(summary.Liver$RDat, annota.Liver$WGT)]
summary.Whole_Blood$gene <- annota.Whole_Blood$ID[match(summary.Whole_Blood$RDat, annota.Whole_Blood$WGT)]

geneEALiver <- unique(intersect(summary.EA$gene, summary.Liver$gene))
geneEAWhole_Blood <- unique(intersect(summary.EA$gene, summary.Whole_Blood$gene))

geneAALiver <- unique(intersect(summary.AA$gene, summary.Liver$gene))
geneAAWhole_Blood <- unique(intersect(summary.AA$gene, summary.Whole_Blood$gene))

## AA Liver
df.hsqAALiver <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneAALiver]), 
                            group=c(rep("P in AA", sum(summary.AA$gene %in% geneAALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneAALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneAALiver)) ),
                            stringsAsFactors = F)

df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))

write.xlsx(df.hsqAALiver, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5b-1", row.names = FALSE, append = TRUE)


## AA Whole_Blood
df.hsqAAWhole_Blood <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneAAWhole_Blood]), 
                                  group=c(rep("P in AA", sum(summary.AA$gene %in% geneAAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
write.xlsx(df.hsqAAWhole_Blood, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5b-2", row.names = FALSE, append = TRUE)



## EA Liver
df.hsqEALiver <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneEALiver]), 
                            group=c(rep("P in EA", sum(summary.EA$gene %in% geneEALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneEALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneEALiver)) ),
                            stringsAsFactors = F)

df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
write.xlsx(df.hsqEALiver, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5b-3", row.names = FALSE, append = TRUE)


## EA Whole_Blood
df.hsqEAWhole_Blood <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneEAWhole_Blood]), 
                                  group=c(rep("P in EA", sum(summary.EA$gene %in% geneEAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
write.xlsx(df.hsqEAWhole_Blood, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig5.xlsx",
           sheetName = "5b-4", row.names = FALSE, append = TRUE)

