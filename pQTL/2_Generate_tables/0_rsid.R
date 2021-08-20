


library(readr)


lookup <- readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")
rsidlookup <- lookup$rsid
names(rsidlookup) <- lookup$SNPid
snpidlookup <- lookup$SNPid
names(snpidlookup) <- lookup$rsid


##################
a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_2.0.txt")
a <- a[,c(1:10,14:15,11:12,16,13)]
a$TopSNP <- rsidlookup[a$TopSNP]
print(a,width=2000)
write_tsv(a,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/1_pQTL_summary_cleaned_2.0_rsid.txt")

a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_2.0.txt")
a <- a[,c(1:10,14:15,11:12,16,13)]
a$TopSNP <- rsidlookup[a$TopSNP]
print(a,width=2000)
write_tsv(a,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/1_pQTL_summary_cleaned_2.0_rsid.txt")

##################
a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis_1.0.txt")
a <- a[,c(1:11,18:19,12:17)]
a$TopSNP <- rsidlookup[a$TopSNP]
a$SNP <- rsidlookup[a$SNP]
print(a,width=2000)
write_tsv(a,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/Tables/2_conditional_analysis_1.0_rsid.txt")

a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis_1.0.txt")
a <- a[,c(1:11,18:19,12:17)]
a$TopSNP <- rsidlookup[a$TopSNP]
a$SNP <- rsidlookup[a$SNP]
print(a,width=2000)
write_tsv(a,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/2_conditional_analysis_1.0_rsid.txt")

##################
a <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/4_replicated_AASK_2.0.txt")
a$replicated <- a$fdr<0.05
a <- a[,c(1:3,5:12,19,14,15,18)]
colnames(a)
a$TopSNP <- rsidlookup[a$TopSNP]
a$AASK_SNP <- rsidlookup[a$AASK_SNP]
print(a,width=2000)
write_tsv(a,"/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/Tables/4_replicated_AASK_2.0_rsid.txt")

