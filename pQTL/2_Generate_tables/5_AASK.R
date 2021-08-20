

## results from revision_500Kb/4_AASK

rm(list=ls())

library(readr)
library(bigreadr)
library(stringr)

### summary
library(readr)
res <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/AASK/replication/results.txt")
dim(res)
res1 <- res[(res$AASK_R2_with_ARIC_topSNP>0.8) & (res$AASK_R2_with_ARIC_topSNP!=1),]; dim(res1)
res2 <- res[res$AASK_R2_with_ARIC_topSNP==1,]; dim(res2)
res3 <- res[res$AASK_R2_with_ARIC_topSNP<=0.8,]; dim(res3)

res <- rbind(res2,res1)
res$fdr <- p.adjust(res$p, method = "fdr")
mean(res$fdr < 0.05) # 0.6905099
mean(res$fdr < 0.01) # 0.6133144
res3$fdr <- NA

mean(p.adjust(res2$p, method = "fdr") < 0.05)

tab1 <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned.txt")

tab12 <- tab1[match(res$SeqID,tab1$SOMAmer),]
mean(sign(tab12$Beta)==sign(res$beta)) # 0.9270538
cor(tab12$Beta,res$beta) # 0.929498
## replicated ones
cor(tab12$Beta[res$fdr<0.05],res$beta[res$fdr<0.05]) # 0.9701473
mean(sign(tab12$Beta[res$fdr<0.05])==sign(res$beta[res$fdr<0.05])) # 0.9969231

tab13 <- tab1[match(res3$SeqID,tab1$SOMAmer),]

tab4 <- rbind(cbind(tab12[,c("SOMAmer","uniprot_id","entrezgenesymbol","Chr","TopSNP","Beta","Pval")], res[,c("AASK_SNP","AASK_R2_with_ARIC_topSNP","beta","p","fdr")]),
              cbind(tab13[,c("SOMAmer","uniprot_id","entrezgenesymbol","Chr","TopSNP","Beta","Pval")], res3[,c("AASK_SNP","AASK_R2_with_ARIC_topSNP","beta","p","fdr")]))

ALT <- character()
AF <- numeric()
i=0
for (chr in 1:22){
  info <- read_tsv(paste0("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/vcf/Black/info/chr",chr,".info.txt"), col_types = cols())
  maf <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/MAF/chr",chr,".afreq"), col_types = cols())
  m <- which(tab4$Chr == chr)
  for (j in 1:length(m)){
    i <- i+1; print(i)

    ALT[i] <- info$ALT[info$ID == tab4$TopSNP[m[j]]]
    tmp <- maf[maf$ID == tab4$TopSNP[m[j]],]
    if(tmp$ALT == ALT[i]){
      AF[i] <- tmp$ALT_FREQS
    }else{
      AF[i] <- 1-tmp$ALT_FREQS
    }
  }
}

tab4$A1 <- ALT
tab4$A1_AF <- AF
tab4$Same_directionaligy <- sign(tab4$Beta)==sign(tab4$beta)

tab4 <- tab4[,c("SOMAmer","uniprot_id","entrezgenesymbol","Chr","TopSNP","A1","A1_AF","Beta","Pval","beta","p","fdr","Same_directionaligy","AASK_SNP","AASK_R2_with_ARIC_topSNP")]

write_tsv(tab4,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/4_replicated_AASK_1.0.txt")


library(readr)

dat <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/4_replicated_AASK_1.0.txt")

dat$sign1=ifelse(sign(dat$Beta)==1, "+","-")
dat$sign2=ifelse(sign(dat$beta)==1, "+","-")
dat$samedir <- paste0(dat$sign1,"/",dat$sign2)
write_tsv(dat, "/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/4_replicated_AASK_2.0.txt")

dat <- dat[1:1400,]
mean(dat$sign1==dat$sign2) # 0.9278571
