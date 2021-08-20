
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



length(RDat.white) #1350

length(RDat.black) #1394

length(intersect(RDat.white, RDat.black)) #1109

load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/genelist/genelist.rds")
length(intersect(paste0(pGene.w,".wgt.RDat"), RDat.white))
length(intersect(paste0(pGene.b,".wgt.RDat"), RDat.black))


