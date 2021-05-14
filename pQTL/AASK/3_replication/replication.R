
rm(list=ls())

library(readr)
library(plink2R)
library(bigreadr)
library(stringr)

lookup <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")
rsidlookup <- lookup$rsid
names(rsidlookup) <- lookup$SNPid
snpidlookup <- lookup$SNPid
names(snpidlookup) <- lookup$rsid

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt')

tab1 <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")
tab1 <- tab1[tab1$SOMAmer %in% seqid,]

seqid <- tab1$SOMAmer
topSNP <- tab1$TopSNP
chr <- tab1$Chr

std <- character(length = length(seqid))
AASK_SNP <- character(length = length(seqid))

R2_withtopSNP <- numeric(length = length(seqid))
Estimate <- numeric(length = length(seqid))
Std.Error <- numeric(length = length(seqid))
P.value <- numeric(length = length(seqid))

for (i in 1:length(seqid)){

  geno <- read_plink(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/window2M/byseq/",seqid[i]), impute="avg")

  if(topSNP[i] %in% geno$bim$V2){
    R2_withtopSNP[i] <- 1 # !!
    std[i] <- "ARIC top SNP is in AASK."
    aasksnp <- topSNP[i]
  }else{
    snp <- c(topSNP[i], geno$bim$V2)
    snprsid <- rsidlookup[snp]
    m <- is.na(snprsid)
    snprsid <- snprsid[!m]
    writeLines(snprsid, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i],".snp"))

    if(m[1]){
      std[i] <- "ARIC top SNP is not in AASK. ARIC top SNP do not have rsid!"
      print(std[i])
      next
    }

    system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink",
                  " --bfile /dcl01/chatterj/data/jzhang2/1000G/GRCh38/AFR/chr",chr[i],
                  " --extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i],".snp",
                  " --make-bed",
                  " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i]),
           ignore.stdout = TRUE)

    system(paste0("/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink1/plink",
                  " --bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i],
                  " --r bin4 --threads 1",
                  " --out /dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i]),
           ignore.stdout = TRUE)

    snp_plink <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i],".bim"),stringsAsFactors = F)$V2

    LD <- readBin(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/tmp/",seqid[i],".ld.bin"),
                  what="numeric", size=4, n=(length(snp_plink))^2)
    LD <- matrix(LD, nrow = length(snp_plink))

    m <- which(snp_plink==snprsid[1])
    if(length(m)==0){
      std[i] <- "ARIC top SNP is not in AASK. ARIC top SNP does not in 1000G!"
      print(std[i])
      next
    }
    R2 <- LD[m,]; R2 <- R2^2
    tmp <- order(R2,decreasing = T)
    tmp <- tmp[-which(tmp==m)]
    m <- tmp[1]
    R2_withtopSNP[i] <- R2[m] # !!
    aasksnp <- snpidlookup[snp_plink[m]]
    names(aasksnp) <- NULL

    std[i] <- paste0("ARIC top SNP is not in AASK. AASK SNPs\' max R2=",R2_withtopSNP[i],".")
  }
  AASK_SNP[i] <- aasksnp #!!

  x <- geno$bed[,aasksnp]
  pheno <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum/invrankpheno/50/", seqid[i], ".pheno"))
  pheno$ID <- paste0(pheno$V1,":",pheno$V2)

  id <- intersect(pheno$ID, names(x))
  dat <- cbind(y=pheno$V3[match(id,pheno$ID)], x=x[match(id,names(x))])
  dat <- as.data.frame(dat)
  fit <- lm(y~.,dat)
  fit <- summary(fit)$coefficients["x",]
  Estimate[i] <- fit["Estimate"]
  Std.Error[i] <-  fit["Std. Error"]
  P.value[i] <-  fit["Pr(>|t|)"]

  print(paste0(i,": ", std[i], " p-value: ", P.value[i]))

}

save(seqid, std, topSNP, AASK_SNP, R2_withtopSNP, Estimate, Std.Error, P.value,
     file=paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/results.RData"))

res <- data.frame(SeqID = seqid,
                  ARIC_topSNP = topSNP, AASK_SNP=AASK_SNP,
                  AASK_R2_with_ARIC_topSNP = R2_withtopSNP,
                  beta= Estimate, p=P.value)
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/replication/results.txt"))
