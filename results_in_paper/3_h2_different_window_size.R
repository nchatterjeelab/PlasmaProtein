
rm(list=ls())
library(stringr)

load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/heritability/h2.Rdat")

sum(res_1M$p_1M<0.01)
sum(res_2M$p_2M<0.01)
a <- intersect(res_1M$SeqID[res_1M$p_1M<0.01], res_2M$SeqID[res_2M$p_2M<0.01])
length(a)
h1 <- res_1M$hsq_1M[res_1M$SeqID %in% a]
h2 <- res_2M$hsq_2M[res_2M$SeqID %in% a]
summary(h1)
summary(h2)


sum(res_1M$p_1M<0.05)
sum(res_2M$p_2M<0.05)
a <- intersect(res_1M$SeqID[res_1M$p_1M<0.05], res_2M$SeqID[res_2M$p_2M<0.05])
length(a)
h1 <- res_1M$hsq_1M[res_1M$SeqID %in% a]
h2 <- res_2M$hsq_2M[res_2M$SeqID %in% a]
summary(h1)
summary(h2)



load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/para1/invrank/heritability/h2.Rdat")

sum(res_1M$p_1M<0.01)
sum(res_2M$p_2M<0.01)
a <- intersect(res_1M$SeqID[res_1M$p_1M<0.01], res_2M$SeqID[res_2M$p_2M<0.01])
length(a)
h1 <- res_1M$hsq_1M[res_1M$SeqID %in% a]
h2 <- res_2M$hsq_2M[res_2M$SeqID %in% a]
summary(h1)
summary(h2)

sum(res_1M$p_1M<0.05)
sum(res_2M$p_2M<0.05)
a <- intersect(res_1M$SeqID[res_1M$p_1M<0.05], res_2M$SeqID[res_2M$p_2M<0.05])
length(a)
h1 <- res_1M$hsq_1M[res_1M$SeqID %in% a]
h2 <- res_2M$hsq_2M[res_2M$SeqID %in% a]
summary(h1)
summary(h2)

a <- readr::read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/Plasma_Protein.pos")
setdiff(res_1M$SeqID[res_1M$p_1M<0.01], paste0(unlist(strsplit(a$WGT, ".wgt.RDat")),".hsq")) # SeqId_17355_56.hsq

#tmp <- c("INHBB","ITIH1","SPP1","BTN3A3","C11orf68","B3GAT3","INHBC","SNUPN")
#tmp <- c("IL1RN","SPP1","BTN3A3","INHBC")
#res_1M[res_1M$SeqID %in% paste0(unlist(strsplit(a$WGT[a$ID %in% tmp], ".wgt.RDat")), ".hsq"), ]

