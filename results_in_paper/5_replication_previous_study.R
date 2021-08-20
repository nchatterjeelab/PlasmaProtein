

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





##################################################
## Replication in ARIC


rm(list=ls())

library(readr)
library(dplyr)
library(plyr)

#lookup <- readRDS("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_nodup.rds")
#rsidlookup <- lookup$rsid
#names(rsidlookup) <- lookup$SNPid
#snpidlookup <- lookup$SNPid
#names(snpidlookup) <- lookup$rsid

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/previous_studies.RData")
dat1 <- dat1[dat1$SOMAmerID %in% seqid,] # 476
dat1 <- dat1[order(as.integer(dat1$CHR)),]
dat2 <- dat2[dat2$SOMAmerID %in% seqid,] # 315
dat2 <- dat2[order(as.integer(dat2$CHR)),]
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
load("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/previous_studies.RData")
dat1 <- dat1[dat1$SOMAmerID %in% seqid,] # 476
dat1 <- dat1[order(as.integer(dat1$CHR)),]
dat2 <- dat2[dat2$SOMAmerID %in% seqid,] # 315
dat2 <- dat2[order(as.integer(dat2$CHR)),]
colnames(dat1) <- paste0("AGES-RS_", colnames(dat1))
colnames(dat2) <- paste0("INTERVAL_", colnames(dat2))
a <- intersect(dat1$'AGES-RS_SOMAmerID', dat2$'INTERVAL_SOMAmerID')
dat <- rbind.fill(cbind(dat1[match(a,dat1$'AGES-RS_SOMAmerID'),], dat2[match(a,dat2$'INTERVAL_SOMAmerID'),]),
                  dat1[!(dat1$'AGES-RS_SOMAmerID' %in% a), ], dat2[!(dat2$'INTERVAL_SOMAmerID' %in% a), ])

seqid <- dat$INTERVAL_SOMAmerID
seqid[is.na(seqid)] <- dat$'AGES-RS_SOMAmerID'[is.na(seqid)]
topSNP1 <- dat$'AGES-RS_TopSNP'
topSNP2 <- dat$INTERVAL_TopSNP
chr <- dat$INTERVAL_CHR
chr[is.na(chr)] <- dat$'AGES-RS_CHR'[is.na(chr)]

std1 <- character(length = length(seqid))
ARIC_SNP1 <- character(length = length(seqid))
R2_withtopSNP1 <- numeric(length = length(seqid))
testedA1 <- character(length = length(seqid))
Estimate1 <- numeric(length = length(seqid))
Std.Error1 <- numeric(length = length(seqid))
P.value1 <- character(length = length(seqid))

std2 <- character(length = length(seqid))
ARIC_SNP2 <- character(length = length(seqid))
R2_withtopSNP2 <- numeric(length = length(seqid))
testedA2 <- character(length = length(seqid))
Estimate2 <- numeric(length = length(seqid))
Std.Error2 <- numeric(length = length(seqid))
P.value2 <- character(length = length(seqid))

for (i in 1:length(seqid)){
  
  sumdata_ARIC <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",seqid[i],".PHENO1.glm.linear"), col_types = "iiccccdciddddc")
  
  ###########################
  ### topSNP1 (AGES-RS)
  if(is.na(topSNP1[i])){
    R2_withtopSNP1[i] <- NA
    Estimate1[i] <- NA
    Std.Error1[i] <- NA
  }else{
    
    if(topSNP1[i] %in% sumdata_ARIC$ID){
      R2_withtopSNP1[i] <- 1 # !!
      std1[i] <- "AGES-RS top SNP is in ARIC."
      aricsnp1 <- topSNP1[i]
      
      tmp <- sumdata_ARIC[sumdata_ARIC$ID==aricsnp1,]
      if(tmp$A1 != dat$'AGES-RS_A1'[i]){ testedA1[i] <- dat$'AGES-RS_A1'[i]; Estimate1[i] <- -tmp$BETA }else{ testedA1[i] <- dat$'AGES-RS_A1'[i]; Estimate1[i] <- tmp$BETA }
      Std.Error1[i] <-  tmp$SE
      P.value1[i] <-  tmp$P
      
    }else{ #print(i)}}
      
      snprsid <- c(topSNP1[i], sumdata_ARIC$ID)
      m <- is.na(snprsid)
      snprsid <- snprsid[!m]
      writeLines(snprsid, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".snp"))
      
      system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2",
                    " --bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh38/AFR/chr",chr[i],
                    " --extract /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".snp",
                    " --make-bed",
                    " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i]),
             ignore.stdout = TRUE)
      
      snp_plink <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".bim"),stringsAsFactors = F)$V2
      m <- which(snp_plink==snprsid[1])
      if(length(m)==0){
        std1[i] <- "AGES-RS top SNP is not in ARIC and 1000G ref."
        print(std1[i])
        R2_withtopSNP1[i] <- NA
        Estimate1[i] <- NA
        Std.Error1[i] <- NA
        
      }else{
        
        system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink",
                      " --bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],
                      " --r bin4 --threads 1",
                      " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i]),
               ignore.stdout = TRUE)
        
        LD <- readBin(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".ld.bin"),
                      what="numeric", size=4, n=(length(snp_plink))^2)
        LD <- matrix(LD, nrow = length(snp_plink))
        
        
        R2 <- LD[m,]; R2 <- R2^2
        tmp <- order(R2,decreasing = T)
        tmp <- tmp[-which(tmp==m)]
        m <- tmp[1]
        R2_withtopSNP1[i] <- R2[m] # !!
        aricsnp1 <- snp_plink[m]
        
        std1[i] <- paste0("AGES-RS top SNP is not in ARIC. ARIC SNPs\' max R2=",R2_withtopSNP1[i],".")
        
        tmp <- sumdata_ARIC[sumdata_ARIC$ID==aricsnp1,]
        
        testedA1[i] <- tmp$A1
        Estimate1[i] <- tmp$BETA
        Std.Error1[i] <-  tmp$SE
        P.value1[i] <-  tmp$P
      }
    }
    ARIC_SNP1[i] <- aricsnp1 #!!
    
  }
  
  ###########################
  ### topSNP2 (INTEVAL)
  if(is.na(topSNP2[i])){
    R2_withtopSNP2[i] <- NA
    Estimate2[i] <- NA
    Std.Error2[i] <- NA
  }else{
    if(topSNP2[i] %in% sumdata_ARIC$ID){
      R2_withtopSNP2[i] <- 1 # !!
      std2[i] <- "INTERVAL top SNP is in ARIC."
      aricsnp2 <- topSNP2[i]
      
      tmp <- sumdata_ARIC[sumdata_ARIC$ID==aricsnp2,]
      if(tmp$A1 != dat$INTERVAL_A1[i]){ testedA2[i] <- dat$INTERVAL_A1[i]; Estimate2[i] <- -tmp$BETA }else{ testedA2[i] <- dat$INTERVAL_A1[i]; Estimate2[i] <- tmp$BETA }
      Std.Error2[i] <-  tmp$SE
      P.value2[i] <-  tmp$P
      
    }else{ #print(i)}}
      
      snprsid <- c(topSNP2[i], sumdata_ARIC$ID)
      m <- is.na(snprsid)
      snprsid <- snprsid[!m]
      writeLines(snprsid, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".snp"))
      
      system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink",
                    " --bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh38/AFR/chr",chr[i],
                    " --extract /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".snp",
                    " --make-bed",
                    " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i]),
             ignore.stdout = TRUE)
      
      snp_plink <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".bim"),stringsAsFactors = F)$V2
      m <- which(snp_plink==snprsid[1])
      if(length(m)==0){
        std2[i] <- "INTERVAL top SNP is not in ARIC and 1000G ref."
        print(std2[i])
        R2_withtopSNP2[i] <- NA
        Estimate2[i] <- NA
        Std.Error2[i] <- NA
      }else{
        
        system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink",
                      " --bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],
                      " --r bin4 --threads 1",
                      " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i]),
               ignore.stdout = TRUE)
        
        LD <- readBin(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".ld.bin"),
                      what="numeric", size=4, n=(length(snp_plink))^2)
        LD <- matrix(LD, nrow = length(snp_plink))
        
        R2 <- LD[m,]; R2 <- R2^2
        tmp <- order(R2,decreasing = T)
        tmp <- tmp[-which(tmp==m)]
        m <- tmp[1]
        R2_withtopSNP2[i] <- R2[m] # !!
        aricsnp2 <- snp_plink[m]
        
        std2[i] <- paste0("INTERVAL top SNP is not in ARIC. ARIC SNPs\' max R2=",R2_withtopSNP2[i],".")
        
        tmp <- sumdata_ARIC[sumdata_ARIC$ID==aricsnp2,]
        
        testedA2[i] <- dat$INTERVAL_A1[i]
        Estimate2[i] <- tmp$BETA
        Std.Error2[i] <-  tmp$SE
        P.value2[i] <-  tmp$P
      }
    }
  }
  ARIC_SNP2[i] <- aricsnp2 #!!
  
  print(paste0(i,": ", P.value1[i], "; ", P.value2[i]))
  
}

save(seqid, dat,
     std1, topSNP1, ARIC_SNP1, R2_withtopSNP1, testedA1, Estimate1, Std.Error1, P.value1,
     std2, topSNP2, ARIC_SNP2, R2_withtopSNP2, testedA2, Estimate2, Std.Error2, P.value2,
     file=paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results.RData"))

res <- data.frame(dat[,1:9],
                  topSNP1 = topSNP1, ARIC_SNP1=ARIC_SNP1, R2_with_topSNP1=R2_withtopSNP1, testedA1_ARIC=testedA1, beta1=Estimate1, se1=Std.Error1, p1=P.value1,
                  dat[,10:18],
                  topSNP2 = topSNP2, ARIC_SNP2=ARIC_SNP2, R2_with_topSNP2=R2_withtopSNP2, testedA2_ARIC=testedA2, beta2=Estimate2, se2=Std.Error1, p2=P.value2)
write.table(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results.txt"),
            quote=F, sep = "\t", row.names = F)




################################################
## summary table
# library(readr)
# replication <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results.txt", col_types = "ccicccddcccdcddccciccccdcccdcddc")
# m <- ( as.numeric(replication$p1) < 1.92*10^(-11)) | ( as.numeric(replication$p2) < 1.5*10^(-11))
# replication <- replication[!is.na(m),]
# m <- m[!is.na(m)]
# mean(m) # 0.988189

# annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
# seqid <- replication$INTERVAL_SOMAmerID
# seqid[is.na(seqid)] <- replication$AGES.RS_SOMAmerID[is.na(seqid)]
# replication$SOMAmer <- seqid
# replication$geneID <- annota$entrezgenesymbol[match(replication$SOMAmer, annota$seqid_in_sample)]
# replication$repl1 <- as.numeric(replication$p1) < 4.92*10^(-5)
# replication$repl2 <- as.numeric(replication$p2) < 4.92*10^(-5)
# tmp <- replication[[2]]; replication[[2]][is.na(tmp)] <- replication[[18]][is.na(tmp)]
# replication <- replication[,c(33,2,34,4,7,9,11,12,14,16,35,20,23,25,27,28,30,32,36)]
# write_tsv(replication, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/4_replicated_previous_studies_1.0.txt"))

# mean((replication$repl1) | (replication$repl2)) # 0.992126
# sum((replication$repl1) | (replication$repl2)) # 504


## summary table
library(readr)
replication <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/4_replicated_previous_studies_1.0.txt")
replication$repl1 <- as.numeric(replication$p1) < 6.89*10^(-5)
replication$repl2 <- as.numeric(replication$p2) < 6.89*10^(-5)
write_tsv(replication, paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/4_replicated_previous_studies_2.0.txt"))

mean((replication$repl1) | (replication$repl2)) # 0.992126
sum((replication$repl1) | (replication$repl2)) # 504




