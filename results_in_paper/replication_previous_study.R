
rm(list=ls())

library(readr)
library(dplyr)
library(plyr)

lookup <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup_nodup.rds")
rsidlookup <- lookup$rsid
names(rsidlookup) <- lookup$SNPid
#snpidlookup <- lookup$SNPid
#names(snpidlookup) <- lookup$rsid

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/previous_studies.RData")
dat1 <- dat1[dat1$SOMAmerID %in% seqid,] # 581
dat1 <- dat1[order(as.integer(dat1$CHR)),]
dat2 <- dat2[dat2$SOMAmerID %in% seqid,] # 504
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

  sumdata_ARIC <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",seqid[i],".PHENO1.glm.linear"), col_types = "iicccccidddcc")
  tmp <- rsidlookup[sumdata_ARIC$ID]; names(tmp) <- NULL
  sumdata_ARIC$rsid <- tmp
  sumdata_ARIC <- sumdata_ARIC[!is.na(tmp),]

  ###########################
  ### topSNP1 (AGES-RS)
  if(is.na(topSNP1[i])){
    R2_withtopSNP1[i] <- NA
    Estimate1[i] <- NA
    Std.Error1[i] <- NA
  }else{


    if(topSNP1[i] %in% sumdata_ARIC$rsid){
      R2_withtopSNP1[i] <- 1 # !!
      std1[i] <- "AGES-RS top SNP is in ARIC."
      aricsnp1 <- topSNP1[i]

      tmp <- sumdata_ARIC[sumdata_ARIC$rsid==aricsnp1,]
      if(tmp$A1 != dat$'AGES-RS_A1'[i]){ testedA1[i] <- dat$'AGES-RS_A1'[i]; Estimate1[i] <- -tmp$BETA }else{ testedA1[i] <- dat$'AGES-RS_A1'[i]; Estimate1[i] <- tmp$BETA }
      Std.Error1[i] <-  tmp$SE
      P.value1[i] <-  tmp$P

    }else{ #print(i)}}

      snprsid <- c(topSNP1[i], sumdata_ARIC$rsid)
      m <- is.na(snprsid)
      snprsid <- snprsid[!m]
      writeLines(snprsid, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".snp"))

      system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink",
                    " --bfile /dcl01/chatterj/data/jzhang2/1000G/GRCh38/AFR/chr",chr[i],
                    " --extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".snp",
                    " --make-bed",
                    " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i]),
             ignore.stdout = TRUE)

      snp_plink <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".bim"),stringsAsFactors = F)$V2
      m <- which(snp_plink==snprsid[1])
      if(length(m)==0){
        std1[i] <- "AGES-RS top SNP is not in ARIC and 1000G ref."
        print(std1[i])
        R2_withtopSNP1[i] <- NA
        Estimate1[i] <- NA
        Std.Error1[i] <- NA

      }else{

        system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink",
                      " --bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],
                      " --r bin4 --threads 1",
                      " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i]),
               ignore.stdout = TRUE)

        LD <- readBin(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/1_",seqid[i],".ld.bin"),
                      what="numeric", size=4, n=(length(snp_plink))^2)
        LD <- matrix(LD, nrow = length(snp_plink))


        R2 <- LD[m,]; R2 <- R2^2
        tmp <- order(R2,decreasing = T)
        tmp <- tmp[-which(tmp==m)]
        m <- tmp[1]
        R2_withtopSNP1[i] <- R2[m] # !!
        aricsnp1 <- snp_plink[m]

        std1[i] <- paste0("AGES-RS top SNP is not in ARIC. ARIC SNPs\' max R2=",R2_withtopSNP1[i],".")

        tmp <- sumdata_ARIC[sumdata_ARIC$rsid==aricsnp1,]

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
    if(topSNP2[i] %in% sumdata_ARIC$rsid){
      R2_withtopSNP2[i] <- 1 # !!
      std2[i] <- "INTERVAL top SNP is in ARIC."
      aricsnp2 <- topSNP2[i]

      tmp <- sumdata_ARIC[sumdata_ARIC$rsid==aricsnp2,]
      if(tmp$A1 != dat$INTERVAL_A1[i]){ testedA2[i] <- dat$INTERVAL_A1[i]; Estimate2[i] <- -tmp$BETA }else{ testedA2[i] <- dat$INTERVAL_A1[i]; Estimate2[i] <- tmp$BETA }
      Std.Error2[i] <-  tmp$SE
      P.value2[i] <-  tmp$P

    }else{ #print(i)}}

      snprsid <- c(topSNP2[i], sumdata_ARIC$rsid)
      m <- is.na(snprsid)
      snprsid <- snprsid[!m]
      writeLines(snprsid, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".snp"))

      system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink",
                    " --bfile /dcl01/chatterj/data/jzhang2/1000G/GRCh38/AFR/chr",chr[i],
                    " --extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".snp",
                    " --make-bed",
                    " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i]),
             ignore.stdout = TRUE)

      snp_plink <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".bim"),stringsAsFactors = F)$V2
      m <- which(snp_plink==snprsid[1])
      if(length(m)==0){
        std2[i] <- "INTERVAL top SNP is not in ARIC and 1000G ref."
        print(std2[i])
        R2_withtopSNP2[i] <- NA
        Estimate2[i] <- NA
        Std.Error2[i] <- NA
      }else{

        system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink",
                      " --bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],
                      " --r bin4 --threads 1",
                      " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i]),
               ignore.stdout = TRUE)

        LD <- readBin(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/tmp/2_",seqid[i],".ld.bin"),
                      what="numeric", size=4, n=(length(snp_plink))^2)
        LD <- matrix(LD, nrow = length(snp_plink))

        R2 <- LD[m,]; R2 <- R2^2
        tmp <- order(R2,decreasing = T)
        tmp <- tmp[-which(tmp==m)]
        m <- tmp[1]
        R2_withtopSNP2[i] <- R2[m] # !!
        aricsnp2 <- snp_plink[m]

        std2[i] <- paste0("INTERVAL top SNP is not in ARIC. ARIC SNPs\' max R2=",R2_withtopSNP2[i],".")

        tmp <- sumdata_ARIC[sumdata_ARIC$rsid==aricsnp2,]

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
     file=paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results.RData"))

res <- data.frame(dat[,1:9],
                  topSNP1 = topSNP1, ARIC_SNP1=ARIC_SNP1, R2_with_topSNP1=R2_withtopSNP1, testedA1_ARIC=testedA1, beta1=Estimate1, se1=Std.Error1, p1=P.value1,
                  dat[,10:19],
                  topSNP2 = topSNP2, ARIC_SNP2=ARIC_SNP2, R2_with_topSNP2=R2_withtopSNP2, testedA2_ARIC=testedA2, beta2=Estimate2, se2=Std.Error1, p2=P.value2)
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results.txt"))



################################################
## summary table
replication <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results.txt", col_types = "ccicccddcccdcddccciccccddcccdcddc")
replication$fdr1 <- p.adjust(as.numeric(replication$p1), method = "fdr")
replication$fdr2 <- p.adjust(as.numeric(replication$p2), method = "fdr")
m <- (replication$fdr1<0.05)|(replication$fdr2<0.05)
replication <- replication[!is.na(m),]
m <- m[!is.na(m)]
mean(m) # 0.9938931

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
seqid <- replication$INTERVAL_SOMAmerID
seqid[is.na(seqid)] <- replication$AGES.RS_SOMAmerID[is.na(seqid)]
replication$SOMAmer <- seqid
replication$geneID <- annota$entrezgenesymbol[match(replication$SOMAmer, annota$seqid_in_sample)]
write_tsv(replication, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Replication_previous_study/results_1.0.txt"))



