

##################################################
## Clean existing studies data

library(readxl)
library(dplyr)
library(stringr)
library(readr)

####################

dat1 <- read_excel("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/pQTLstudies/pQTL_summary_AGES-RS.xlsx",
                   sheet = "Supplementary Table S1")
dat1 <- dat1[dat1$Chr_SNP == dat1$Chr_protein,]
dat1$BP <- as.integer(gsub(",","",dat1$BP))
dat1$`Target Gene_start` <- as.integer(gsub(",","",dat1$`Target Gene_start`))
dat1 <- dat1[abs(dat1$BP - dat1$`Target Gene_start`)<0.5*10^6,]
dat1 <- dat1[!str_detect(dat1$Uniprot, " |,"),]
dat1 <- dat1[as.numeric(dat1$`P-value`) < 1.92*10^(-11),]
dat1 <- dat1[(dat1$A1_freq>0.01) & (dat1$A1_freq<0.99),]
dat1 <- dat1[nchar(dat1$A1)==1,]
dat1 <- dat1[nchar(dat1$A2)==1,]
dat1$SOMAmerID <- paste0("SeqId_", gsub("-","_",unlist(lapply(str_split(dat1$`Apatamer 2`,"_"),FUN=function(x){paste0(x[1],"_",x[2])}))))

SOMAmerID <- unique(dat1$SOMAmerID)
res <- tibble()
for(i in 1:length(SOMAmerID)){
  tmp <- dat1[dat1$SOMAmerID==SOMAmerID[i],]
  res <- rbind(res,tmp[which.min(tmp$`P-value`),])
}
dat1 <- res
dat1 <- dat1[,c("SOMAmerID","Uniprot","Chr_protein","rsID","BP","A1","BETA","SE","P-value")]
colnames(dat1) <- c("SOMAmerID","UniProtID","CHR","TopSNP","POS","A1","BETA","SE","P")

####################

annota <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt")
a <- paste0("chr",annota$chromosome_name,":",annota$transcription_start_site,"-",annota$transcription_start_site)
#writeLines(a, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/tmp/dat2_hg19.bed")
b1 <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/tmp/dat2_hg38.bed")
b2 <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/tmp/dat2_hg38.err")
annota <- annota[-which(a %in% b2),]
annota$TSS_hg19 <- as.integer(unlist(lapply(str_split(b1, ":|-"), FUN=function(x){x[2]})))

dat2_soma <- dir("/dcs04/nilanjan/data/jzhang2/INTERVAL_proteome/Somalogic")
dat2_soma <- data.frame(indx=dat2_soma, seqid=paste0("SeqId_", unlist(lapply(str_split(dat2_soma,"\\."),FUN=function(x){paste0(x[2],"_",x[3])}))))
dat2_soma <- inner_join(dat2_soma, annota[,c(1,8:10)], by=c("seqid"="seqid_in_sample"))
colnames(dat2_soma)[3:4] <- c("CHR","TSS_hg38")
dat2_soma <- dat2_soma[order(dat2_soma$CHR),]


#for (i in 1:nrow(dat2_soma)){
#  for (chr in 1:22){
#    tmp <- paste0("/",dat2_soma$indx[i],"_chrom_",chr,"_meta_final_v1.tsv.gz")
#    if( !(paste0(dat2_soma$indx[i],"_chrom_",chr,"_meta_final_v1.tsv.gz") %in% dir( paste0("/dcs04/nilanjan/data/jzhang2/INTERVAL_proteome/Somalogic/",dat2_soma$indx[i],"/") )) ){
#      print(paste0(dat2_soma$indx[i],"  ","chr",chr))
#    }
#  }
#}

progBar <- function(ii, N, per = 10) {
  #ii is current iteration.
  #N is total number of iterations to perform.
  #per is step (percent) when to update the progress bar. We need this as when multiple iterations are being performed progress bar gets updated too often and output is messed up.
  if (ii %in% seq(1, N, per)) {
    x <- round(ii * 100 / N)
    message("[ ",
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""),
            " ] ", x, "%", "\r",
            appendLF = FALSE)
    if (ii == N) cat("\r")
  }
}

library(doMC)
library(foreach)
registerDoMC(12)

res <- tibble()
for (chr in 1:22){

  print(chr)

  bim <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh37/EUR/chr",chr,".bim"))
  maf <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh37/EUR/chr",chr,"_freq.afreq"))
  dat2_soma_tmp <- dat2_soma[dat2_soma$CHR==chr,]

  niter <- nrow(dat2_soma_tmp)
  tmp <- foreach(i=1:niter, ii = icount(), .combine='rbind') %dopar% {
    tmp <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/INTERVAL_proteome/Somalogic/",dat2_soma_tmp$indx[i],"/",dat2_soma_tmp$indx[i],"_chrom_",chr,"_meta_final_v1.tsv.gz"))
    tmp <- tmp[abs(tmp$position - dat2_soma_tmp$TSS_hg19[i])<0.5*10^6,]
    tmp <- tmp[nchar(tmp$Allele1)==1,]
    tmp <- tmp[nchar(tmp$Allele2)==1,]
    tmp <- inner_join(tmp,bim[,c(2,4)],by=c("position"="V4"))
    tmp <- inner_join(tmp,maf[,c(2,4,5)],by=c("V2"="ID"))
    tmp <- tmp[(tmp$ALT_FREQS>0.01) & (tmp$ALT_FREQS<0.99),]
    progBar(ii, niter, per=5)
    tmp$SOMAmerID <- dat2_soma_tmp$seqid[i]
    tmp <- tmp[which.min(tmp$`log(P)`),]
    return(tmp)
  }
  res <- rbind(res, tmp)
}
save(res,file="/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/dat2_res.RData")

load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/dat2_res.RData")
res1 <- res[,c(2,9,3,4,6,7)]
colnames(res1) <- c("CHR","TopSNP","POS","A1","BETA","SE")
annota <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt")
dat2 <- cbind(SOMAmerID=res$SOMAmerID, #? res$SOMAmerID
              UniProtID=annota$uniprot_id[match(res$SOMAmerID, annota$seqid_in_sample)],
              res1,
              P=exp(as.numeric(res[,"log(P)"])), stringsAsFactors=F)
dat2 <- dat2[dat2$P<1.5*10^(-11),]
rownames(dat2) <- NULL
dat2$A1 <- toupper(dat2$A1)

save(dat1, dat2, file="/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/previous_studies.RData")

