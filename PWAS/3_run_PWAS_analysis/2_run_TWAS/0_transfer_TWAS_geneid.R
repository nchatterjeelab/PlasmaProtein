


#
#tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")
#ensembl <- character()
#for (tissue in tissue_list) {
#  ensembl <- c(ensembl,suppressMessages(read.table(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR/GTEXv8.",tissue,".pos"), header = T, stringsAsFactors = F))$ID)
#  ensembl <- unique(ensembl)
#  print(tissue)
#}
#for (tissue in tissue_list) {
#  ensembl <- c(ensembl,suppressMessages(read.table(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue,".pos"), header = T, stringsAsFactors = F))$ID)
#  ensembl <- unique(ensembl)
#  print(tissue)
#}
#writeLines(ensembl, "/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/ensemblid.txt")

#library(biomaRt)
#library(stringr)
#library(readr)
#ensemblid <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/ensemblid.txt")
#ensemblid <- unlist(lapply(str_split(ensemblid,"[.]"), FUN = function (x){x[1]}))
#ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#se.dat <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
#                filters = "ensembl_gene_id",
#                values = ensemblid,
#                mart = ensembl)
#se.dat <- se.dat[se.dat$hgnc_symbol!="",]
#dim(se.dat)
#write_tsv(se.dat, "/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/convertgene/dictionary.txt")
#a <- setdiff(ensemblid, se.dat$ensembl_gene_id)
#writeLines(a, "/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/convertgene/nogeneid.txt")
#
#library("org.Hs.eg.db")
#library(clusterProfiler)
##bitr("ABO", fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
#b <- bitr(a, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
#b <- bitr(a, fromType="ENSEMBL", toType=, OrgDb="org.Hs.eg.db")



library(dplyr)
library(readr)
library(stringr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")
#se.dat <- read.table("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/gencode.v26.GRCh38.genes.gtf", stringsAsFactors = F, sep = "\t")
#tmp <- str_split(se.dat$V9, "; ")
#genename <- unlist(lapply(tmp, function (x){ x[str_detect(x,"gene_name")] })); genename <- gsub("gene_name ","",genename)
#emsemblname <- unlist(lapply(tmp, function (x){ x[str_detect(x,"gene_id")] })); emsemblname <- gsub("gene_id ","",emsemblname)
#se.dat <- data.frame(emsemblname,genename,stringsAsFactors = F)
#se.dat <- unique(se.dat)
#write_tsv(se.dat, "/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/convertgene/gtexv8dictionary.txt")
se.dat <- read_tsv("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/convertgene/gtexv8dictionary.txt")


for (j in 1:length(tissue_list)) {
  tissue <- tissue_list[j]
  annot.t <- suppressMessages(read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue,".pos")))
  #annot.t <- suppressMessages(read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR/GTEXv8.",tissue,".pos")))
  annot.t <- annot.t[!is.na(annot.t$CHR),]
  gene.t <- annot.t$ID
  #gene.t <- unlist(lapply(str_split(gene.t,"[.]"), FUN = function (x){x[1]}))

  m <- match(gene.t,se.dat$emsemblname)
  annot.t$ID <- se.dat$genename[m]
  m <- is.na(annot.t$ID)
  annot.t <- annot.t[!m,]
  write_tsv(annot.t, paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue,".geneid.pos"))
  #write_tsv(annot.t, paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR/GTEXv8.",tissue,".geneid.pos"))

}


GTEx_V8_n_gene.ALL <- integer()
for (j in 1:length(tissue_list)) {
  tissue <- tissue_list[j]
  annot.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue,".geneid.pos"))
  GTEx_V8_n_gene.ALL[j] <- nrow(annot.t)
}
names(GTEx_V8_n_gene.ALL) <- tissue_list
saveRDS(GTEx_V8_n_gene.ALL, "/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_n_gene.rds")

