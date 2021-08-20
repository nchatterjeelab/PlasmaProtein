

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

library(readr)
library(plink2R)

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/byseq/"))

if(ethnic == "White"){
  n_peer <- 90
}else{
  n_peer <- 80
}

tab1 <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/1_pQTL_summary_cleaned_2.0_rsid.txt"))

seqid <- tab1$SOMAmer

flagPAV <- rep(0, length(seqid))
flagnoPAV <- rep(0, length(seqid))

Estimate <- numeric()
Std.Error <- numeric()
P.value <- numeric()
for (i in 1:length(seqid)){
  coding <- read.table(paste0("/dcs04/nilanjan/data/diptavo/coding_info/",ethnic,"/",seqid[i],".txt"), header = T, stringsAsFactors = F)
  PAV <- coding$SNP[coding$coding==1]
  PAV <- unique(PAV)

  if(length(PAV)==0){ flagnoPAV[i] <- 2; print("noPAV"); next }
  writeLines(PAV, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/byseq/",seqid[i],".PAV"))

  if(length(PAV)==1){
    PAV.prune.in <- PAV
  }else{
    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink",
                  " --keep-allele-order",
                  " --bfile /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/window1M/byseq/",seqid[i],
                  " --extract /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/byseq/",seqid[i],".PAV",
                  " --threads 1",
                  " --indep-pairwise 100 10 0.9",
                  " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/byseq/",seqid[i]),
    ignore.stdout = TRUE, ignore.stderr = TRUE)

    PAV.prune.in <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/byseq/",seqid[i],".prune.in"))
    tmp <- table(PAV.prune.in)
    PAV.prune.in <- names(tmp[tmp==1])
  }

  topSNP <- tab1$TopSNP[i]

  geno <- read_plink(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/window1M/byseq/",seqid[i]), impute="avg")
  dat <- cbind(x=geno$bed[,topSNP], geno$bed[,PAV.prune.in])

  R2 <- cor(dat)^2
  maxR2 <- max(R2[1,-1])

  if(maxR2>0.9){ flagPAV[i] <- 1; print("PAV"); next } else if ( maxR2<0.1 ){ flagnoPAV[i] <- 1; print("noPAV"); next }

  pheno <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum/invrankpheno/",n_peer,"/", seqid[i], ".pheno"))
  pheno$ID <- paste0(pheno$V1,":",pheno$V2)

  id <- intersect(pheno$ID, rownames(dat))
  dat <- cbind(y=pheno$V3[match(id,pheno$ID)], dat[match(id,rownames(dat)),])
  dat <- as.data.frame(dat)
  fit <- lm(y~.,dat)
  fit <- summary(fit)$coefficients["x",]
  Estimate[i] <- fit["Estimate"]
  Std.Error[i] <-  fit["Std. Error"]
  P.value[i] <-  fit["Pr(>|t|)"]
  print(i)
}
save(seqid, flagPAV, flagnoPAV, Estimate, Std.Error, P.value,
     file=paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/results.RData"))


thres <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/peernum_permutation/",n_peer,"/permutation/permutations_all.thresholds.txt"))
res <- data.frame(SeqID = seqid,
                  flagPAV = flagPAV,flagnoPAV=flagnoPAV,
                  thres = thres$V2[match(seqid,thres$V1)],
                  beta= Estimate, p=P.value,
                  marg_beta=tab1$Beta, margp=tab1$Pval)
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/results.txt"))


### Clean results
res <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/epitope_binding_effect/results.txt"))
tmp <- rep("0.1 <= R2 <= 0.9", nrow(res))
tmp[as.logical(res$flagPAV)] <- "R2 > 0.9"
tmp[as.logical(res$flagnoPAV)] <- "R2 < 0.1"
res$PAVLD <- tmp
res$sig <- res$p < res$thres
res$sig[res$PAVLD=="R2 < 0.1"] <- TRUE
write_tsv(res,paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic,"/pQTL/Tables/3_PAV.txt"))
