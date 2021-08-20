
################################################################
################################################################
################################################################

rm(list=ls())

library(stringr)
library(readr)

bimid <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M_pre/byseq")
bimid <- bimid[str_detect(bimid,".bim")]

res <- character()
res1 <- character()
for (i in 1:length(bimid)){
  tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M_pre/byseq/", bimid[i]),
                    col.names = c('Chr', 'ID', 'bpm', 'Pos_b38', 'A1', 'A2'),
                    stringsAsFactors=F)
  res1 <- c(res1, tmp$ID)
  tmp0 <- format(tmp$Pos_b38, trim = T, scientific = F)
  tmp <- paste0("chr", tmp$Chr, ":", tmp0, "-", tmp0)
  res <- c(res, tmp)
  print(i)
}
df <- data.frame(SNP=res1, GRCh38.key=res, stringsAsFactors = F)
df <- unique(df)
write_tsv(df, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_GRCh38_ID.txt")
writeLines(unique(res), "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_GRCh38_cleaned.txt")

#
#GRCh38 <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_GRCh38_cleaned.txt")
#mapped <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_GRCh37_cleaned.bed")
#err <- readLinGes("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_GRCh37_cleaned_err.bed")
#GRCh38 <- GRCh38[!(GRCh38 %in% err)]
#
#lookup <- data.frame(stringsAsFactors = F,
#                     Chr = str_sub(gsub(":.*", "", mapped), start = 4),
#                     Pos_b38 = gsub(".*-", "", GRCh38),
#                     Pos_b37 = gsub(".*-", "", mapped),
#                     GRCh38.key = GRCh38,
#                     GRCh37.key = mapped)
#write_tsv(lookup, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_lookup_cleaned.txt")
#
#library(dplyr)
#lookup <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_lookup_cleaned.txt")
#lookup <- inner_join(lookup,df, by="GRCh38.key")
#write_tsv(lookup, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/SNPconvert/ARIC_lookup_cleaned_ID.txt")
