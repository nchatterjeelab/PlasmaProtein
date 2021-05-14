########################################################################
########################################################################
## PWAS model and h2
########################################################################
########################################################################


rm(list=ls())
library(stringr)

tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
hsq <- tmp[str_detect(tmp,"hsq")]
RDat <- tmp[str_detect(tmp,"RDat")]
hsq.all <- numeric()
rsq.enet <- numeric()
rsq.top1 <- numeric()
pval.hsq <- numeric()
pval.enet <- numeric()
pval.top1 <- numeric()
for (i in 1:length(RDat)){
  load(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/coefs_remove_ambiguous_snp/", RDat[i]))
  hsq.all[i] <- hsq[1]
  pval.hsq[i] <- hsq.pv
  rsq.enet[i] <- cv.performance[1,"enet"]
  rsq.top1[i] <- cv.performance[1,"top1"]
  pval.enet[i] <- cv.performance[2,"enet"]
  pval.top1[i] <- cv.performance[2,"top1"]
  print(i)
}
summary.black <- list(M=length(RDat),
                      RDat = unlist(lapply(str_split(RDat, '[.]wgt[.]RDat'), FUN=function (x){x[1]})),
                      hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
saveRDS(summary.black, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/summary.black.rds")
RDat.black <- RDat

tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
hsq <- tmp[str_detect(tmp,"hsq")]
RDat <- tmp[str_detect(tmp,"RDat")]
hsq.all <- numeric()
rsq.enet <- numeric()
rsq.top1 <- numeric()
pval.hsq <- numeric()
pval.enet <- numeric()
pval.top1 <- numeric()
for (i in 1:length(RDat)){
  load(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp/", RDat[i]))
  hsq.all[i] <- hsq[1]
  pval.hsq[i] <- hsq.pv
  rsq.enet[i] <- cv.performance[1,1]
  rsq.top1[i] <- cv.performance[1,2]
  pval.enet[i] <- cv.performance[2,1]
  pval.top1[i] <- cv.performance[2,2]
  print(i)
}
summary.white <-  list(M=length(RDat),
                       RDat = unlist(lapply(str_split(RDat, '[.]wgt[.]RDat'), FUN=function (x){x[1]})),
                       hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
saveRDS(summary.white, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/summary.white.rds")
RDat.white <- RDat


####################################
### summary table

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/1_Heritability/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/1_Heritability/summary.white.rds")


df.AA <- data.frame(seqid=summary.AA$RDat, 
                    h2=summary.AA$hsq.all,
                    R2=summary.AA$rsq.enet, 
                    stringsAsFactors = F)
df.EA <- data.frame(seqid=summary.EA$RDat, 
                    h2=summary.EA$hsq.all, 
                    R2=summary.EA$rsq.enet, 
                    stringsAsFactors = F)
res <- full_join(df.EA, df.AA, by="seqid")

overlapping <- intersect(summary.AA$RDat, summary.EA$RDat)
rest1 <- setdiff(summary.EA$RDat,overlapping)
rest2 <- setdiff(summary.AA$RDat,overlapping)

res <- res[match(c(overlapping,rest1,rest2),res$seqid),]
colnames(res) <- c("seqid","EA_h2","EA_R2","AA_h2","AA_R2")

annota <- read_tsv("/Users/jnz/Document/JHU/Research/ARIC/Abbreviated annotation visits 3 and 5.txt")

res <- left_join(res, annota[,c("seqid_in_sample","uniprot_id","entrezgenesymbol","target","targetfullname")],
                   by=c("seqid"="seqid_in_sample"))
write_tsv(res,"/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Tables/protein_h2.txt")

########################################################################
########################################################################
## TWAS model and h2
########################################################################
########################################################################

rm(list=ls())
library(stringr)

tmp <- dir("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS")
tmp <- tmp[str_detect(tmp,"hsq")]
tissue_list <- gsub("[.].*","",tmp)

for(t in 1:length(tissue_list)){

  RDat <- dir(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/",tissue_list[t]))
  RDat <- RDat[str_detect(RDat,"RDat")]
  hsq.all <- numeric()
  rsq.enet <- numeric()
  rsq.top1 <- numeric()
  pval.hsq <- numeric()
  pval.enet <- numeric()
  pval.top1 <- numeric()
  for (i in 1:length(RDat)){
    load(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/",tissue_list[t],"/", RDat[i]))
    hsq.all[i] <- hsq[1]
    pval.hsq[i] <- hsq.pv
    rsq.enet[i] <- cv.performance[1,"enet"]
    rsq.top1[i] <- cv.performance[1,"top1"]
    pval.enet[i] <- cv.performance[2,"enet"]
    pval.top1[i] <- cv.performance[2,"top1"]
    print(paste0(t,"_",i))
  }
  summary <- list(RDat=RDat, hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
  saveRDS(summary, paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/RESULTS/summary.",tissue_list[t],".rds"))

}

####################################
###  summary table

rm(list=ls())
library(stringr)

tmp <- dir("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS")
tmp <- tmp[str_detect(tmp,"hsq")]
tissue_list <- gsub("[.].*","",tmp)

med.h2 <- numeric()
med.r2.top1 <- numeric()
med.r2.enet <- numeric()
med.adj.top1 <- numeric()
med.adj.enet <- numeric()
for(t in 1:length(tissue_list)){
  summary <- readRDS(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/RESULTS/summary.",tissue_list[t],".rds"))

  med.h2[t] <- median(summary$hsq.all)
  med.r2.top1[t] <- median(summary$rsq.top1)
  med.r2.enet[t] <- median(summary$rsq.enet)
  med.adj.top1[t] <- median(summary$rsq.top1/summary$hsq.all)
  med.adj.enet[t] <- median(summary$rsq.enet/summary$hsq.all)
  print(t)
}
data.frame(tissue_list,
           round(med.h2,4),
           round(med.r2.top1,4),
           round(med.r2.enet,4),
           round(med.adj.top1,4),
           round(med.adj.enet,4), stringsAsFactors = F)


for(t in 1:length(tissue_list)){
  summary <- readRDS(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/RESULTS/summary.",tissue_list[t],".rds"))

  mean.h2[t] <- mean(summary$hsq.all)
  mean.r2.top1[t] <- mean(summary$rsq.top1)
  mean.r2.enet[t] <- mean(summary$rsq.enet)
  mean.adj.top1[t] <- mean(summary$rsq.top1/summary$hsq.all)
  mean.adj.enet[t] <- mean(summary$rsq.enet/summary$hsq.all)
  print(t)
}
data.frame(tissue_list,
           round(mean.h2,4),
           round(mean.r2.top1,4),
           round(mean.r2.enet,4),
           round(mean.adj.top1,4),
           round(mean.adj.enet,4), stringsAsFactors = F)





