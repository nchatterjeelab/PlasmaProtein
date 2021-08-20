########################################################################
########################################################################
## PWAS model and h2
########################################################################
########################################################################



rm(list=ls())
library(stringr)

tmp <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/coefs")
hsq <- tmp[str_detect(tmp,"hsq")]
RDat <- tmp[str_detect(tmp,"RDat")]
hsq.all <- numeric()
rsq.enet <- numeric()
rsq.top1 <- numeric()
pval.hsq <- numeric()
pval.enet <- numeric()
pval.top1 <- numeric()
for (i in 1:length(RDat)){
  load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/coefs/", RDat[i]))
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
saveRDS(summary.black, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/summary.black.rds")

RDat.black <- RDat

tmp <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs")
hsq <- tmp[str_detect(tmp,"hsq")]
RDat <- tmp[str_detect(tmp,"RDat")]
hsq.all <- numeric()
rsq.enet <- numeric()
rsq.top1 <- numeric()
pval.hsq <- numeric()
pval.enet <- numeric()
pval.top1 <- numeric()
for (i in 1:length(RDat)){
  load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs/", RDat[i]))
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
saveRDS(summary.white, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/summary.white.rds")

RDat.white <- RDat


####################################
### summary table

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")


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
res <- res[,c(1,6:9,2:5)]
res$sigboth <- (!is.na(res$EA_h2)) & (!is.na(res$AA_h2))
write_tsv(res,"/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/5_protein_h2.txt")


## model summary

tmp <- rbind(c(round(median(summary.AA$hsq.all),4),
           round(median(summary.AA$rsq.top1),4),
           round(median(summary.AA$rsq.enet),4),
           round(median(summary.AA$rsq.top1/summary.AA$hsq.all),4),
           round(median(summary.AA$rsq.enet/summary.AA$hsq.all),4)),
      c(round(median(summary.EA$hsq.all),4),
                 round(median(summary.EA$rsq.top1),4),
                 round(median(summary.EA$rsq.enet),4),
                 round(median(summary.EA$rsq.top1/summary.EA$hsq.all),4),
                 round(median(summary.EA$rsq.enet/summary.EA$hsq.all),4)))
View(tmp)

tmp <- rbind(c(round(mean(summary.AA$hsq.all),4),
               round(mean(summary.AA$rsq.top1),4),
               round(mean(summary.AA$rsq.enet),4),
               round(mean(summary.AA$rsq.top1/summary.AA$hsq.all),4),
               round(mean(summary.AA$rsq.enet/summary.AA$hsq.all),4)),
             c(round(mean(summary.EA$hsq.all),4),
               round(mean(summary.EA$rsq.top1),4),
               round(mean(summary.EA$rsq.enet),4),
               round(mean(summary.EA$rsq.top1/summary.EA$hsq.all),4),
               round(mean(summary.EA$rsq.enet/summary.EA$hsq.all),4)))
View(tmp)

## average improvement

round(median((summary.EA$rsq.enet-summary.EA$rsq.top1)/summary.EA$rsq.top1),4) # 0.3563
round(median((summary.AA$rsq.enet-summary.AA$rsq.top1)/summary.AA$rsq.top1),4) # 0.3982

########################################################################
########################################################################
## TWAS model and h2
########################################################################
########################################################################

rm(list=ls())
library(stringr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

for(t in 1:length(tissue_list)){

  RDat <- dir(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue_list[t]))
  RDat <- RDat[str_detect(RDat,"RDat")]
  hsq.all <- numeric()
  rsq.enet <- numeric()
  rsq.top1 <- numeric()
  pval.hsq <- numeric()
  pval.enet <- numeric()
  pval.top1 <- numeric()
  for (i in 1:length(RDat)){
    load(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue_list[t],"/", RDat[i]))
    hsq.all[i] <- hsq[1]
    pval.hsq[i] <- hsq.pv
    rsq.enet[i] <- cv.performance[1,"enet"]
    rsq.top1[i] <- cv.performance[1,"top1"]
    pval.enet[i] <- cv.performance[2,"enet"]
    pval.top1[i] <- cv.performance[2,"top1"]
    print(paste0(t,"_",i))
  }
  summary <- list(RDat=RDat, hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
  saveRDS(summary, paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v7/summary.",tissue_list[t],".rds"))

}


rm(list=ls())
library(stringr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")

for(t in 1:length(tissue_list)){
  
  RDat <- dir(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue_list[t]))
  RDat <- RDat[str_detect(RDat,"RDat")]
  hsq.all <- numeric()
  rsq.enet <- numeric()
  rsq.top1 <- numeric()
  pval.hsq <- numeric()
  pval.enet <- numeric()
  pval.top1 <- numeric()
  for (i in 1:length(RDat)){
    load(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue_list[t],"/", RDat[i]))
    hsq.all[i] <- hsq[1]
    pval.hsq[i] <- hsq.pv
    rsq.enet[i] <- cv.performance[1,"enet"]
    rsq.top1[i] <- cv.performance[1,"top1"]
    pval.enet[i] <- cv.performance[2,"enet"]
    pval.top1[i] <- cv.performance[2,"top1"]
    print(paste0(t,"_",i))
  }
  summary <- list(RDat=RDat, hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
  saveRDS(summary, paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v8.ALL/summary.",tissue_list[t],".rds"))
  
}

####################################
###  summary table

## V7

rm(list=ls())
library(stringr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

med.h2 <- numeric()
med.r2.top1 <- numeric()
med.r2.enet <- numeric()
med.adj.top1 <- numeric()
med.adj.enet <- numeric()
for(t in 1:length(tissue_list)){
  summary <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v7/summary.",tissue_list[t],".rds"))

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

mean.h2 <- numeric()
mean.r2.top1 <- numeric()
mean.r2.enet <- numeric()
mean.adj.top1 <- numeric()
mean.adj.enet <- numeric()
for(t in 1:length(tissue_list)){
  summary <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v7/summary.",tissue_list[t],".rds"))

  mean.h2[t] <- mean(summary$hsq.all)
  mean.r2.top1[t] <- mean(summary$rsq.top1)
  mean.r2.enet[t] <- mean(summary$rsq.enet)
  mean.adj.top1[t] <- mean(summary$rsq.top1/summary$hsq.all)
  mean.adj.enet[t] <- mean(summary$rsq.enet/summary$hsq.all)
  print(t)
}
data.frame(#tissue_list,
           round(mean.h2,4),
           round(mean.r2.top1,4),
           round(mean.r2.enet,4),
           round(mean.adj.top1,4),
           round(mean.adj.enet,4), stringsAsFactors = F)


## V8

rm(list=ls())
library(stringr)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")

med.h2 <- numeric()
med.r2.top1 <- numeric()
med.r2.enet <- numeric()
med.adj.top1 <- numeric()
med.adj.enet <- numeric()
for(t in 1:length(tissue_list)){
  summary <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v8.ALL/summary.",tissue_list[t],".rds"))
  
  med.h2[t] <- median(summary$hsq.all, na.rm = T)
  med.r2.top1[t] <- median(summary$rsq.top1, na.rm = T)
  med.r2.enet[t] <- median(summary$rsq.enet, na.rm = T)
  med.adj.top1[t] <- median(summary$rsq.top1/summary$hsq.all, na.rm = T)
  med.adj.enet[t] <- median(summary$rsq.enet/summary$hsq.all, na.rm = T)
  print(t)
}
data.frame(tissue_list,
           round(med.h2,4),
           round(med.r2.top1,4),
           round(med.r2.enet,4),
           round(med.adj.top1,4),
           round(med.adj.enet,4), stringsAsFactors = F)

mean.h2 <- numeric()
mean.r2.top1 <- numeric()
mean.r2.enet <- numeric()
mean.adj.top1 <- numeric()
mean.adj.enet <- numeric()
for(t in 1:length(tissue_list)){
  summary <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v8.ALL/summary.",tissue_list[t],".rds"))
  
  mean.h2[t] <- mean(summary$hsq.all, na.rm = T)
  mean.r2.top1[t] <- mean(summary$rsq.top1, na.rm = T)
  mean.r2.enet[t] <- mean(summary$rsq.enet, na.rm = T)
  mean.adj.top1[t] <- mean(summary$rsq.top1/summary$hsq.all, na.rm = T)
  mean.adj.enet[t] <- mean(summary$rsq.enet/summary$hsq.all, na.rm = T)
  print(t)
}
data.frame(#tissue_list,
  round(mean.h2,4),
  round(mean.r2.top1,4),
  round(mean.r2.enet,4),
  round(mean.adj.top1,4),
  round(mean.adj.enet,4), stringsAsFactors = F)



########################################################################
########################################################################
## TWAS model #samples and #features
########################################################################
########################################################################

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")
NS <- integer()
NF <- integer()
for (t in 1:length(tissue_list)) {
  tmp <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue_list[t],".P01.pos"), col_types = cols())
  NS[t] <- tmp$N[1]
  NF[t] <- nrow(tmp)
}
res <- data.frame(tissue_list,NS,NF)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")
NS <- integer()
NF <- integer()
for (t in 1:length(tissue_list)) {
  tmp <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue_list[t],".geneid.pos"), col_types = cols())
  NS[t] <- tmp$N[1]
  NF[t] <- nrow(tmp)
}
res <- data.frame(tissue_list,NS,NF)



