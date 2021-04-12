
rm(list=ls())
library(stringr)

RDat <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
RDat <- RDat[str_detect(RDat,"RDat")]
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
  rsq.enet[i] <- cv.performance["rsq","enet"]
  rsq.top1[i] <- cv.performance["rsq","top1"]
  pval.enet[i] <- cv.performance["pval","enet"]
  pval.top1[i] <- cv.performance["pval",'top1']
  print(i)
}
summary.black <- list(M=length(RDat),hsq.all=hsq.all,
                      rsq.enet=rsq.enet, rsq.top1=rsq.top1,
                      pval.hsq=pval.hsq, pval.enet=pval.enet,
                      pval.top1=pval.top1,
                      Soma=unlist(lapply(str_split(RDat,'.wgt.RDat'), FUN=function (x){x[1]})))
saveRDS(summary.black, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/summary.black.rds")
RDat.black <- RDat

RDat <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
RDat <- RDat[str_detect(RDat,"RDat")]
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
  rsq.enet[i] <- cv.performance["rsq","enet"]
  rsq.top1[i] <- cv.performance["rsq","top1"]
  pval.enet[i] <- cv.performance["pval","enet"]
  pval.top1[i] <- cv.performance["pval",'top1']
  print(i)
}
summary.white <-  list(M=length(RDat),hsq.all=hsq.all,
                       rsq.enet=rsq.enet, rsq.top1=rsq.top1,
                       pval.hsq=pval.hsq, pval.enet=pval.enet,
                       pval.top1=pval.top1,
                       Soma=unlist(lapply(str_split(RDat,'.wgt.RDat'), FUN=function (x){x[1]})))
saveRDS(summary.white, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/summary.white.rds")
RDat.white <- RDat

#> length(RDat.white)
#[1] 1319
#> length(RDat.black)
#[1] 1377
#> length(intersect(RDat.white, RDat.black))
#[1] 1098


#pQTL <- read.table( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_permutation/120/conditional/allsig.txt"))
#a <- unique(pQTL$V1); length(a)  #[1] 1605
#b <- unlist(lapply(str_split(RDat.white,'.wgt.RDat'), FUN=function (x){x[1]})); length(b)  #[1] 1397
#length(intersect(a,b))


#RDat <- dir("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/Liver")
#
#
#pQTL <- read.table( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_permutation/",n_peer,"/conditional/allsig.txt"))
#a <- unique(pQTL$V1); length(a)  #[1] 1605
#b <- unlist(lapply(str_split(res_1M$SeqID,".hsq"),FUN=function (x){x[1]})); length(b)  #[1] 1397
#length(intersect(a,b))RDat <- RDat[str_detect(RDat,"RDat")]
#hsq.all <- numeric()
#rsq.enet <- numeric()
#rsq.top1 <- numeric()
#pval.hsq <- numeric()
#pval.enet <- numeric()
#pval.top1 <- numeric()
#for (i in 1:length(RDat)){
#  load(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/Liver/", RDat[i]))
#  hsq.all[i] <- hsq[1]
#  pval.hsq[i] <- hsq.pv
#  rsq.enet[i] <- cv.performance["rsq","enet"]
#  rsq.top1[i] <- cv.performance["rsq","top1"]
#  pval.enet[i] <- cv.performance["pval","enet"]
#  pval.top1[i] <- cv.performance["pval",'top1']
#  print(i)
#}
#summary.Liver <-  list(M=length(RDat),hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
#saveRDS(summary.Liver, "/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/RESULTS/summary.Liver.rds")
#
#
#RDat <- dir("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/Whole_Blood")
#RDat <- RDat[str_detect(RDat,"RDat")]
#hsq.all <- numeric()
#rsq.enet <- numeric()
#rsq.top1 <- numeric()
#pval.hsq <- numeric()
#pval.enet <- numeric()
#pval.top1 <- numeric()
#for (i in 1:length(RDat)){
#  load(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/Whole_Blood/", RDat[i]))
#  hsq.all[i] <- hsq[1]
#  pval.hsq[i] <- hsq.pv
#  rsq.enet[i] <- cv.performance["rsq","enet"]
#  rsq.top1[i] <- cv.performance["rsq","top1"]
#  pval.enet[i] <- cv.performance["pval","enet"]
#  pval.top1[i] <- cv.performance["pval",'top1']
#  print(i)
#}
#summary.Whole_Blood <-  list(M=length(RDat),hsq.all=hsq.all, rsq.enet=rsq.enet, rsq.top1=rsq.top1, pval.hsq=pval.hsq, pval.enet=pval.enet, pval.top1=pval.top1)
#saveRDS(summary.Whole_Blood, "/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/RESULTS/summary.Whole_Blood.rds")

summary.white <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/summary.white.rds")
summary.black <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/summary.black.rds")

data.frame(h2 = c(round(summary(summary.black$hsq.all),4)[3],round(summary(summary.white$hsq.all),4)[3]),
r2_top1= c(round(summary(summary.black$rsq.top1),4)[3],round(summary(summary.white$rsq.top1),4)[3]),
r2_enet= c(round(summary(summary.black$rsq.enet),4)[3],round(summary(summary.white$rsq.enet),4)[3]),
acc_top1= c(round(summary(summary.black$rsq.top1/summary.black$hsq.all),4)[3],round(summary(summary.white$rsq.top1/summary.white$hsq.all),4)[3]),
acc_top1= c(round(summary(summary.black$rsq.enet/summary.black$hsq.all),4)[3],round(summary(summary.white$rsq.enet/summary.white$hsq.all),4)[3]))

data.frame(h2 = c(round(summary(summary.black$hsq.all),4)[4],round(summary(summary.white$hsq.all),4)[4]),
r2_top1= c(round(summary(summary.black$rsq.top1),4)[4],round(summary(summary.white$rsq.top1),4)[4]),
r2_enet= c(round(summary(summary.black$rsq.enet),4)[4],round(summary(summary.white$rsq.enet),4)[4]),
acc_top1= c(round(summary(summary.black$rsq.top1/summary.black$hsq.all),4)[4],round(summary(summary.white$rsq.top1/summary.white$hsq.all),4)[4]),
acc_top1= c(round(summary(summary.black$rsq.enet/summary.black$hsq.all),4)[4],round(summary(summary.white$rsq.enet/summary.white$hsq.all),4)[4]))


#### FURTHER ANALYSIS ON MY MAC

#/GRCh38/PWAS/White/1_invrank/2_results/results_TWAS.R