################################
#Generating indiv file for FAST#
# From ThuyVy
# /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered
################################

library(readr)
library(dplyr)
key <- read_csv("/dcs01/arking/ARIC_static/ARIC_Data/gwa_official_idlist_rev121023.csv")
fam <- read.table("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr1.fam",
                  stringsAsFactors = F)
id.old <- fam$V2
id.old.tomatch <- gsub(".*_","",id.old)
tmp <- data.frame(id.old, id.old.tomatch, stringsAsFactors = F)
matched <- inner_join(tmp, key, by = c("id.old.tomatch"="gwasid"))
write_tsv(matched[,c(1,3)], "/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt", col_names=F)


library(readr)
library(dplyr)
key <- read_csv("/dcs01/arking/ARIC_static/ARIC_Data/gwa_official_idlist_rev121023.csv")
fam <- read.table("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr1.fam",
                  stringsAsFactors = F)
id.old <- fam$V2
id.old.tomatch <- gsub(".*_","",id.old)
tmp <- data.frame(id.old, id.old.tomatch, stringsAsFactors = F)
matched <- inner_join(tmp, key, by = c("id.old.tomatch"="gwasid"))
write_tsv(matched[,c(1,3)], "/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt", col_names=F)



