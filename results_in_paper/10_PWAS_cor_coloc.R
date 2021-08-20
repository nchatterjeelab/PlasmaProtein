library(readxl)
library(dplyr)
library(readr)

urateid <- c("SeqId_13676_46","SeqId_7955_195","SeqId_17692_2","SeqId_19622_7","SeqId_6897_38","SeqId_8307_47","SeqId_15686_49","SeqId_17765_3","SeqId_8900_28","SeqId_8403_18")
urategene <- c("INHBB","ITIH1","BTN3A3","INHBA","B3GAT3","C11orf68","INHBC","SNUPN","NEO1","FASN")

dat1 <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/conditional_analysis/Urate_all-cleaned_samegene.out")
dat2 <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/conditional_analysis/Urate_all-cleaned_samegene_v8.out")
corr <- full_join(dat1[,c(2:6,1,13, 7:11, 14:15)], dat2[,c(2,1, 13, 7:11, 14:15)], by=c("PWAS_hit","tissue"))


coloc <- read_excel("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_06_revision2/Suppl_tables_9Aug2021_DD_JZ.xlsx",
                    sheet = "ST8.2-- coloc(PP.H4) with eQTLs", skip=2)
tmp <- coloc[match(urateid,coloc$`SOMAmer ID`),]
rescoloc <- numeric(); resid <- character(); restissue <- character()
ii=0
for (i in 1:length(urateid)) {
  for (j in 3:ncol(tmp)) {
    ii <- ii+1
    resid[ii] <- as.character(tmp[i,2])
    restissue[ii] <- colnames(tmp)[j]
    rescoloc[ii] <- as.numeric(tmp[i,j])
  }
}
coloc <- data.frame(gene=resid, tissue=restissue,pph4=rescoloc, stringsAsFactors = F)


corr$tissue1 <- paste0(corr$PWAS_hit,"-",unlist(lapply(strsplit(corr$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")})))
coloc$tissue1 <- paste0(coloc$gene,"-",unlist(lapply(strsplit(colnames(tmp)[-1:-2], "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")})))

res <- left_join(corr, coloc[,2:4], by="tissue1")
res$tissue.x <- res$tissue.y
res <- res[,-which(colnames(res) %in% c("tissue.y","tissue1","TWAS_hit.x","TWAS_hit.y"))]
res <- res[( (!is.na(res$Corr_of_hits.x)) | (!is.na(res$Corr_of_hits.y)) ) & !is.na(res$tissue.x),]

res <- res[order(res$tissue.x),]
res <- full_join(data.frame(name=urategene), res,by=c("name"="PWAS_hit"))

df <- res[!(is.na(res$TWAS_p.x)), ]
# df <- tibble()
# for (i in 1:length(urategene)) {
#   tmp <- res[res$name == urategene[i],]
#   tmp <- tmp[!is.na(tmp$TWAS_p.x),]
#   tmp <- tmp[which.min(as.numeric(gsub("[*]","",tmp$TWAS_p.x))),]
#   df <- rbind(df, tmp)
# }
# df <- df[,c(1:6,8:11,7,12:13,15:18,14,19:21)]
write_tsv(df, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/8_PWAS_cor_coloc_urate.txt")


library(readxl)
library(dplyr)
library(readr)

goutid <- c("SeqId_5353_89","SeqId_17692_2","SeqId_15686_49")
goutgene <- c("IL1RN","BTN3A3","INHBC")

dat1 <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/conditional_analysis/Gout_all-cleaned_samegene.out")
dat2 <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/conditional_analysis/Gout_all-cleaned_samegene_v8.out")
corr <- full_join(dat1[,c(2:6,1,13, 7:11, 14:15)], dat2[,c(2,1, 13, 7:11, 14:15)], by=c("PWAS_hit","tissue"))


coloc <- read_excel("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_06_revision2/Suppl_tables_9Aug2021_DD_JZ.xlsx",
                    sheet = "ST8.2-- coloc(PP.H4) with eQTLs", skip=2)
tmp <- coloc[match(goutid,coloc$`SOMAmer ID`),]
rescoloc <- numeric(); resid <- character(); restissue <- character()
ii=0
for (i in 1:length(goutid)) {
  for (j in 3:ncol(tmp)) {
    ii <- ii+1
    resid[ii] <- as.character(tmp[i,2])
    restissue[ii] <- colnames(tmp)[j]
    rescoloc[ii] <- as.numeric(tmp[i,j])
  }
}
coloc <- data.frame(gene=resid, tissue=restissue,pph4=rescoloc, stringsAsFactors = F)


corr$tissue1 <- paste0(corr$PWAS_hit,"-",unlist(lapply(strsplit(corr$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")})))
coloc$tissue1 <- paste0(coloc$gene,"-",unlist(lapply(strsplit(colnames(tmp)[-1:-2], "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")})))

res <- left_join(corr, coloc[,2:4], by="tissue1")
res$tissue.x <- res$tissue.y
res <- res[,-which(colnames(res) %in% c("tissue.y","tissue1","TWAS_hit.x","TWAS_hit.y"))]
res <- res[( (!is.na(res$Corr_of_hits.x)) | (!is.na(res$Corr_of_hits.y)) ) & !is.na(res$tissue.x),]

res <- res[order(res$tissue.x),]
res <- full_join(data.frame(name=goutgene), res,by=c("name"="PWAS_hit"))

df <- res[!(is.na(res$TWAS_p.x)), ]
# df <- tibble()
# for (i in 1:length(urategene)) {
#   tmp <- res[res$name == urategene[i],]
#   tmp <- tmp[!is.na(tmp$TWAS_p.x),]
#   tmp <- tmp[which.min(as.numeric(gsub("[*]","",tmp$TWAS_p.x))),]
#   df <- rbind(df, tmp)
# }
# df <- df[,c(1:6,8:11,7,12:13,15:18,14,19:21)]
write_tsv(df, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/8_PWAS_cor_coloc_gout.txt")


