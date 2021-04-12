#######################################
# split white and black

rm(list=ls())

library(readr)
pc.w <- read_csv("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/White/6_Documentation/aric.f3v2.w.PC.9489.csv")
pc.b <- read.table("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/Black/6_Documentation/aric_f3v2_b_10pc.txt",header=T, stringsAsFactors = F)

sample.w <- read_delim("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/samples.txt", col_names=F, delim=" ")
sample.b <- read_delim("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/samples.txt", col_names=F, delim=" ")

pc.w <- pc.w[match(sample.w$X2, pc.w$indiv_id),]
pc.b <- pc.b[match(sample.b$X2, pc.b$ID),]

colnames(pc.w)[1] <- "SampleId"
colnames(pc.b)[1] <- "SampleId"

write_tsv(pc.w, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pc.txt')
write_tsv(pc.b, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pc.txt')
