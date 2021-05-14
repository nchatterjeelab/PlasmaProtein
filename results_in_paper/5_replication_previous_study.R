library(readxl)
library(dplyr)

dat1 <- read_excel("/Users/jnz/Document/JHU/Research/pQTL/1_2020_4782p_5457i_Human serum proteome profoundly overlaps with genetic signatures of disease.xlsx", 
                   sheet = " ST1; Exome array GWAS results")
dat1 <- dat1[dat1$`Chr. pQTLs` == dat1$`Chr. protein affected`,]
dat1$`hg19_coordinates Gene start` <- gsub(",","",dat1$`hg19_coordinates Gene start`)
dat1$`hg19_coordinates Gene start` <- as.integer(dat1$`hg19_coordinates Gene start`)
dat1 <- dat1[abs(dat1$hg19_pos - dat1$`hg19_coordinates Gene start`)<10^6,]
dat1 <- dat1[!str_detect(dat1$UniProt, " |,"),]
dat1$SOMAmerID <- paste0("SeqId_", gsub("-","_",unlist(lapply(str_split(dat1$SOMAmer,"_"),FUN=function(x){x[1]}))))
dat1 <- dat1[(dat1$AGES_A1_AF>0.01) & (dat1$AGES_A1_AF<0.99),]
dat1 <- dat1[nchar(dat1$A1)==1,]

SOMAmerID <- unique(dat1$SOMAmerID)
res <- tibble()
for(i in 1:length(SOMAmerID)){
  tmp <- dat1[dat1$SOMAmerID==SOMAmerID[i],]
  res <- rbind(res,tmp[which.min(tmp$`P-value`),])
}
dat1 <- res
dat1 <- dat1[,c("SOMAmerID","UniProt","Chr. protein affected","dbSNPID","hg19_pos","A1","BETA","STAT","P-value")]
colnames(dat1) <- c("SOMAmerID","UniProtID","CHR","TopSNP","POS","A1","BETA","STAT","P")

# dat2 <- read_excel("/Users/jnz/Document/JHU/Research/pQTL/pQTL_paper/supplementary_41586_2018_175_MOESM4_ESM.xlsx", 
#                    sheet = "ST4 - pQTL summary", skip=4)
dat2 <- read_csv("/Users/jnz/Document/JHU/Research/pQTL/datasource/pQTL_summary_ATLAS.csv", col_types = "ccccccicccccdcccddcddcddciccccc")
dat2 <- dat2[dat2$`cis/ trans` == "cis",]
dat2 <- dat2[!str_detect(dat2$UniProt, ","),]
dat2$SOMAmerID <- paste0("SeqId_", unlist(lapply(str_split(dat2$`SOMAmer ID`,"\\."),FUN=function(x){paste0(x[2],"_",x[3])})))
dat2 <- dat2[(dat2$EAF>0.01) & (dat2$EAF<0.99),]
dat2 <- dat2[(nchar(dat2$`Effect Allele (EA)`)==1) & (nchar(dat2$`Other Allele (OA)`)==1),]

SOMAmerID <- unique(dat2$SOMAmerID)
res <- tibble()
for(i in 1:length(SOMAmerID)){
  tmp <- dat2[dat2$SOMAmerID==SOMAmerID[i],]
  res <- rbind(res,tmp[which.min(tmp$p_2),])
}
dat2 <- res
dat2 <- dat2[,c("SOMAmerID","UniProt","Chr","Sentinel variant*","Pos","Effect Allele (EA)","Other Allele (OA)","beta_2","SE_2","p_2")]
colnames(dat2) <- c("SOMAmerID","UniProtID","CHR","TopSNP","POS","A1","A0","BETA","SE","P")

save(dat1, dat2, file="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*RData/previous_studies.RData")

