
## white

rm(list=ls())

library(stringr)

tmp <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs")
hsqseq <- tmp[str_detect(tmp,"hsq")]

hsq_1M <- numeric()
hsqp_1M <- numeric()
length(hsqseq)
for (i in 1:length(hsqseq)){
  tmp <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs/",hsqseq[i]))
  hsq_1M[i] <- tmp$V2
  hsqp_1M[i] <- tmp$V4
  print(i)
}
res_1M <- data.frame(SeqID = gsub(".hsq","",hsqseq), hsq_1M = hsq_1M, p_1M = hsqp_1M, stringsAsFactors = F)


tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/h2_2M/")
hsqseq <- tmp[str_detect(tmp,"hsq")]

hsq_2M <- numeric()
hsqp_2M <- numeric()
length(hsqseq)
for (i in 1:length(hsqseq)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/h2_2M/",hsqseq[i]))
  hsq_2M[i] <- tmp$V2
  hsqp_2M[i] <- tmp$V4
  print(i)
}
res_2M <- data.frame(SeqID = gsub(".hsq","",hsqseq), hsq_2M = hsq_2M, p_2M = hsqp_2M, stringsAsFactors = F)

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/heritability/")
save(res_1M, res_2M, file = "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/heritability/h2.Rdat")

res_1M <- res_1M[(res_1M$p_1M < 0.01) & (res_1M$hsq_1M >=0),]
res_2M <- res_2M[(res_2M$p_2M < 0.01) & (res_2M$hsq_2M >=0),]
dim(res_1M)
#[1] 1350    3
dim(res_2M)
#[1] 1461    3
length(intersect(res_1M$SeqID,res_2M$SeqID))
#[1] 1262

library(dplyr)
res <- inner_join(res_1M, res_2M, by="SeqID")
summary(res$hsq_1M)
summary(res$hsq_2M)
library(MASS)
wilcox.test(res$hsq_1M, res$hsq_2M, paired=TRUE)



#########################
#########################

## black

rm(list=ls())

library(stringr)

tmp <- dir("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/coefs")
hsqseq <- tmp[str_detect(tmp,"hsq")]

hsq_1M <- numeric()
hsqp_1M <- numeric()
length(hsqseq)
for (i in 1:length(hsqseq)){
  tmp <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/coefs/",hsqseq[i]))
  hsq_1M[i] <- tmp$V2
  hsqp_1M[i] <- tmp$V4
  print(i)
}
res_1M <- data.frame(SeqID = gsub(".hsq","",hsqseq), hsq_1M = hsq_1M, p_1M = hsqp_1M, stringsAsFactors = F)


tmp <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/h2_2M/")
hsqseq <- tmp[str_detect(tmp,"hsq")]

hsq_2M <- numeric()
hsqp_2M <- numeric()
length(hsqseq)
for (i in 1:length(hsqseq)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/h2_2M/",hsqseq[i]))
  hsq_2M[i] <- tmp$V2
  hsqp_2M[i] <- tmp$V4
  print(i)
}
res_2M <- data.frame(SeqID = gsub(".hsq","",hsqseq), hsq_2M = hsq_2M, p_2M = hsqp_2M, stringsAsFactors = F)

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/heritability/")
save(res_1M, res_2M, file = "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/heritability/h2.Rdat")

res_1M <- res_1M[(res_1M$p_1M < 0.01) & (res_1M$hsq_1M >=0),]
res_2M <- res_2M[(res_2M$p_2M < 0.01) & (res_2M$hsq_2M >=0),]
dim(res_1M)
#[1] 1350    3
dim(res_2M)
#[1] 1461    3
length(intersect(res_1M$SeqID,res_2M$SeqID))
#[1] 1262

library(dplyr)
res <- inner_join(res_1M, res_2M, by="SeqID")
summary(res$hsq_1M)
summary(res$hsq_2M)
library(MASS)
wilcox.test(res$hsq_1M, res$hsq_2M, paired=TRUE)



