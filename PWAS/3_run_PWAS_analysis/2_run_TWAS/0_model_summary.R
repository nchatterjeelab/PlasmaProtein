
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



