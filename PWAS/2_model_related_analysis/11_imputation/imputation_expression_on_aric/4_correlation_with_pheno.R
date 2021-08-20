
rm(list = ls())

suppressMessages(library(readr))

################################################################################
## true phenotype (inverse-rank normalized)

pheno_ID <- read.table("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/imputation_expression_on_aric/geno/chr1.fam")$V2

prot <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos")
pheno_matrix = matrix(nrow=7213,ncol=nrow(prot))
seqid <- unlist(lapply(prot$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)}))
for (i in 1:nrow(prot)) {
  pheno <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/90/",seqid[i],".pheno"))
  pheno <- pheno[match(pheno_ID,pheno$V2),]
  pheno_matrix[,i] = pheno$V3
  print(i)
}
colnames(pheno_matrix) <- seqid
write_tsv(as.data.frame(pheno_matrix), paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/imputation_expression_on_aric/aric_pheno.txt"))



################################################################################
## correlation between imputed gene expression and true phenotype (inverse-rank normalized)

library(dplyr)
library(readr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

annot.p <- suppressMessages(read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos"))
imputed.p <- suppressMessages(read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/imputation_expression_on_aric/aric_pheno.txt")))
m <- match(unlist(lapply(annot.p$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)})),
           colnames(imputed.p))
imputed.p <- imputed.p[,m]

gene.p <- annot.p$ID

res <- tibble()
for (tissue in tissue_list) {
  #if ( tissue %in% c("Brain_Hypothalamus","Skin_Sun_Exposed_Lower_leg","Whole_Blood") )
  #  next

  annot.t <- suppressMessages(read_tsv(paste0('/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/',tissue,'.P01.pos')))
  imputed.t <- suppressMessages(read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/imputation_expression_on_aric/imputed/",tissue,".txt")))
  gene.t <- annot.t$ID

  gene <- intersect(unique(gene.p), unique(gene.t))
  cor <- numeric()
  for (i in 1:length(gene)){
    imputed.p_tmp <- imputed.p[,annot.p$ID == gene[i]]
    imputed.t_tmp <- imputed.t[,annot.t$ID == gene[i]]

    if(ncol(imputed.p_tmp)>1){
      tmp <- imputed.p_tmp[1]
      for (j in 2:ncol(imputed.p_tmp)) {
        tmp <- tmp+imputed.p_tmp[j]
      }
      imputed.p_tmp <- tmp
    }

    if(ncol(imputed.t_tmp)>1){
      tmp <- imputed.t_tmp[1]
      for (j in 2:ncol(imputed.t_tmp)) {
        tmp <- tmp+imputed.t_tmp[j]
      }
      imputed.t_tmp <- tmp
    }

    cor[i] <- as.numeric(cor(imputed.p_tmp, imputed.t_tmp))

  }

  res <- rbind(res,
               data.frame(tissue=tissue, cor= cor, gene=gene))

  rm(list = c("annot.t", "imputed.t","gene.t"))

  print(tissue)

}

write_tsv(res,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/imputation_expression_on_aric/correlations.txt")




res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/correlations.txt")
summ <- res %>% group_by(tissue) %>%
    summarise(ave = mean(cor, na.rm=T),
              med = median(cor, na.rm = T),
              N=length(gene))
summ <- summ[order(summ$med,decreasing=T),]

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE",
                          sep = '\t', comment.char = '', stringsAsFactors = F)
a <- unlist(lapply(strsplit(summ$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
summ$tissue <- gtex.colors$V1[match(a, b)]


