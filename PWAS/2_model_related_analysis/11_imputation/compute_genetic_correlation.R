##########################################

## V7

library(dplyr)
library(readr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

annot.p <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos"), col_types = cols())
imputed.p <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_prot.txt"), col_types = cols())
m <- match(unlist(lapply(annot.p$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)})),
           colnames(imputed.p))
imputed.p <- imputed.p[,m]

gene.p <- annot.p$ID

res <- tibble()
for (tissue in tissue_list) {
  annot.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue,".P01.pos"), col_types = cols())
  imputed.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v7/",tissue,".txt"), col_types = cols())
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
  res <- rbind(res, data.frame(tissue=tissue, cor= cor, gene=gene))
  
  print(tissue)
  
}
write_tsv(res,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/correlations_v7.txt")


##########################################

## V8

library(dplyr)
library(readr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")

annot.p <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos"), col_types = cols())
imputed.p <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_prot.txt"), col_types = cols())
m <- match(unlist(lapply(annot.p$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)})),
           colnames(imputed.p))
imputed.p <- imputed.p[,m]

gene.p <- annot.p$ID

res <- tibble()
for (tissue in tissue_list) {
  annot.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue,".geneid.pos"), col_types = cols())
  imputed.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v8.ALL/",tissue,".txt"), col_types = cols())
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
  res <- rbind(res, data.frame(tissue=tissue, cor= cor, gene=gene))

  print(tissue)

}
write_tsv(res,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/correlations_v8.ALL.txt")

##########################################
 ## clean results

res <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/correlations_v7.txt")
summ <- res %>% group_by(tissue) %>%
    summarise(ave = mean(cor, na.rm=T),
              med = median(cor, na.rm = T),
              N=length(gene))
summ <- summ[order(summ$med,decreasing=T),]

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")
NS <- integer()
NF <- integer()
for (t in 1:length(tissue_list)) {
  tmp <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue_list[t],".P01.pos"), col_types = cols())
  NS[t] <- tmp$N[1]
  NF[t] <- nrow(tmp)
}
Nsample_v7 <- data.frame(tissue_list,NS,NF)

summ$Nsample <- Nsample_v7$NS[match(summ$tissue,Nsample_v7$tissue_list)]

summ_v7 <- summ
colnames(summ_v7) <- paste0("v7_",colnames(summ_v7))


res <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/correlations_v8.ALL.txt")
summ <- res %>% group_by(tissue) %>%
    summarise(ave = mean(cor, na.rm=T),
              med = median(cor, na.rm = T),
              N=length(gene))

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")
NS <- integer()
NF <- integer()
for (t in 1:length(tissue_list)) {
  tmp <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue_list[t],".geneid.pos"), col_types = cols())
  NS[t] <- tmp$N[1]
  NF[t] <- nrow(tmp)
}
Nsample_v8 <- data.frame(tissue_list,NS,NF)

summ$Nsample <- Nsample_v8$NS[match(summ$tissue,Nsample_v8$tissue_list)]

summ_v8 <- summ
colnames(summ_v8) <- paste0("v8_",colnames(summ_v8))


newtissue <- setdiff(summ_v8$v8_tissue, summ_v7$v7_tissue)
tmp1 <- summ_v8[!(summ_v8$v8_tissue %in% newtissue), ]; tmp2 <- summ_v8[summ_v8$v8_tissue %in% newtissue, ]
summ <- plyr::rbind.fill(cbind(summ_v7, tmp1[match(summ_v7$v7_tissue, tmp1$v8_tissue),]),
                         tmp2[order(tmp2$v8_med,decreasing=T),])
summ <- summ[,c(6,2,3,5,4,7,8,10,9)]
#colnames(summ)[1] <- "tissue"

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE",
                          sep = '\t', comment.char = '', stringsAsFactors = F)
a <- unlist(lapply(strsplit(summ$v8_tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
summ$tissue <- gtex.colors$V1[match(a, b)]
if(sum(is.na(summ$tissue))>0){ summ$tissue[is.na(summ$tissue)] <- a[is.na(summ$tissue)]}
write_tsv(summ,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/supp_tab.txt")




##########################################

## imputed expression on ARIC

library(dplyr)
library(readr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

annot.p <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos"), col_types = cols())
imputed.p <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_prot.txt"), col_types = cols())
m <- match(unlist(lapply(annot.p$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)})),
           colnames(imputed.p))
imputed.p <- imputed.p[,m]

gene.p <- annot.p$ID

res <- tibble()
for (tissue in tissue_list) {
  annot.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue,".P01.pos"), col_types = cols())
  imputed.t <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v7/",tissue,".txt"), col_types = cols())
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
  res <- rbind(res, data.frame(tissue=tissue, cor= cor, gene=gene))

  print(tissue)

}
write_tsv(res,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/correlations_v7.txt")



#tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")
#a <- dir(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v7/"))
#for (tissue in tissue_list) {
#  if( !(paste0(tissue,".txt") %in% a) ){
#    print(tissue)
#  }
#}
