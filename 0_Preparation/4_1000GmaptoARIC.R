
################################################################
################################################################
################################################################

rm(list=ls())
library(readr)
library(stringr)
library(dplyr)

ethnic <- "EUR"
ethnic1 <- "White"

ARIC <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_ID.txt"))
a <- str_split(ARIC$GRCh38.key,"chr|:|-")
ARIC$Chr <- as.integer(unlist(lapply(a,FUN=function (x){x[2]})))
ARIC$Pos <- unlist(lapply(a,FUN=function (x){x[3]}))
write_tsv(ARIC, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_map.txt"))
# 6181856
ARIC <- read_tsv( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_map.txt"))
RSID_lookup <- tibble()
for (chr in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/1000G/GRCh38/",ethnic,"/chr", chr,".bim"),
                    col.names = c('Chr', 'Rsid', 'bpm', 'Pos_b38', 'A1', 'A2'),
                    stringsAsFactors=F)
  tmp1 <- ARIC[ARIC$Chr == chr,]
  tmp1 <- tmp1[!(tmp1$Pos %in% as.integer(names(table(tmp1$Pos)[table(tmp1$Pos)>1]))),]
  tmp <- tmp[tmp$Pos_b38 %in% tmp1$Pos, ]
  tmp <- tmp[!(tmp$Pos_b38 %in% as.integer(names(table(tmp$Pos_b38)[table(tmp$Pos_b38)>1]))),]
  tmp1 <- inner_join(tmp1, tmp[,c(2,4)], by=c("Pos"="Pos_b38"))
  RSID_lookup <- rbind(RSID_lookup, tmp1)
  print(chr)
}
# 6074650 rows
writeLines(RSID_lookup$Rsid, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_RSID.txt"))
writeLines(RSID_lookup$SNP, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_ARICID.txt"))
write_tsv(RSID_lookup, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_RSID_lookup.txt"))


rm(list=ls())
library(readr)
library(stringr)
library(dplyr)

ethnic <- "AFR"
ethnic1 <- "Black"

ARIC <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_ID.txt"))
a <- str_split(ARIC$GRCh38.key,"chr|:|-")
ARIC$Chr <- as.integer(unlist(lapply(a,FUN=function (x){x[2]})))
ARIC$Pos <- unlist(lapply(a,FUN=function (x){x[3]}))
write_tsv(ARIC, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_map.txt"))
# 8228168
ARIC <- read_tsv( paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_map.txt"))
RSID_lookup <- tibble()
for (chr in 1:22){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/1000G/GRCh38/",ethnic,"/chr", chr,".bim"),
                    col.names = c('Chr', 'Rsid', 'bpm', 'Pos_b38', 'A1', 'A2'),
                    stringsAsFactors=F)
  tmp1 <- ARIC[ARIC$Chr == chr,]
  tmp1 <- tmp1[!(tmp1$Pos %in% as.integer(names(table(tmp1$Pos)[table(tmp1$Pos)>1]))),]
  tmp <- tmp[tmp$Pos_b38 %in% tmp1$Pos, ]
  tmp <- tmp[!(tmp$Pos_b38 %in% as.integer(names(table(tmp$Pos_b38)[table(tmp$Pos_b38)>1]))),]
  tmp1 <- inner_join(tmp1, tmp[,c(2,4)], by=c("Pos"="Pos_b38"))
  RSID_lookup <- rbind(RSID_lookup, tmp1)
  print(chr)
} # 8070219 rows
writeLines(RSID_lookup$Rsid, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_RSID.txt"))
writeLines(RSID_lookup$SNP, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_ARICID.txt"))
write_tsv(RSID_lookup, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/2M_ARIC_GRCh38_RSID_lookup.txt"))


################################################################
################################################################
################################################################


rm(list=ls())

ethnic <- "AFR"
ethnic1 <- "Black"

dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/LDref/",ethnic))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/3_1000GmaptoARIC/",ethnic))

for (chr in 1:22){

b <- paste0("#!/usr/bin/env bash
#$ -N LDref_", chr,"
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G
#$ -m e

/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--bfile /dcl01/chatterj/data/jzhang2/1000G/GRCh38/",ethnic,"/chr",chr," \\
--extract /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/ARIC_GRCh38_RSID.txt \\
--make-bed \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/LDref/",ethnic,"/chr",chr,"


")
  print(chr)

  writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/0_Preparation/3_1000GmaptoARIC/",ethnic,"/chr", chr,".sh"))

}



################################################################
################################################################
################################################################

rm(list=ls())
library(readr)
library(stringr)
library(dplyr)

ethnic <- "AFR"
ethnic1 <- "Black"

dat <- read_tsv(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/ARIC_GRCh38_RSID_lookup.txt"))
write_tsv(dat[,c("SNP","Rsid")], paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/",ethnic1,"/PWAS/SNPconvert/ARIC_to_RSID.txt"),col_names = F)
