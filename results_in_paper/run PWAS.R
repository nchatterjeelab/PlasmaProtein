########################################################################

# Urate (urate summary data has RSID column)

system("cp /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Urate_EA_cleaned_rsid_BETA_SE_TWAS.txt /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Urate_EA_cleaned_rsid_BETA_SE_PWAS.txt")

# Gout (gout summary data doesn't have RSID column. Matching SNP ID is needed. Performed on MARCC.)

folder='CKDGen'
summname='gout_chr1_22_LQ_IQ06_mac10_EA_171.tbl.nstud9.summac400'
respath1='Gout_EA'

Allele1='Allele1'
Allele2='Allele2'
Effect='Effect'
StdErr='StdErr'
type='tsv'
Pval='P-value'
SNPID='MarkerName'

summ <- read_tsv(paste0('/home-2/jzhan218@jhu.edu/work/jzhan218/summary_data/', folder,'/', summname))


summ1 <- summ

a <- str_split(summ$MarkerName,"_")

summ$Chr <- as.integer(unlist(lapply(a,FUN=function (x){x[1]})))
summ$Pos <- as.integer(unlist(lapply(a,FUN=function (x){x[2]})))

summ <- summ[,c("Chr", "Pos", Allele1, Allele2, Effect, StdErr)]
colnames(summ) <- c('Chr', 'Pos', 'A1', 'A2', 'BETA', 'SE')
summ$BETA <- as.numeric(summ$BETA)
summ$SE <- as.numeric(summ$SE)
summ <- summ[(!is.na(summ$BETA)) & (!is.na(summ$SE)),]
summ$Z <- summ$BETA/summ$SE

lookup <- read_tsv("/home-2/jzhan218@jhu.edu/work/jzhan218/PWAS/Results/White/ARIC_lookup_cleaned_RSID.txt")

tmp0 <- format(summ$Pos, trim = T, scientific = F)
summ$GRCh37.key <- paste0("chr", summ$Chr, ":", tmp0, "-", tmp0)
summ <- inner_join(summ, lookup[c("SNP","GRCh37.key")], by="GRCh37.key")
summ <- summ[,c('SNP', 'A1', 'A2', 'Z','BETA','SE')]
write_tsv(summ, paste0('/home-2/jzhan218@jhu.edu/work/jzhan218/summary_data/', folder,'/', respath1, '_cleaned_rsid_BETA_SE.txt'))



########################################################################
## clean table

library(dplyr)
library(readr)

results <- tibble()
for (chr in 1:22) {
    results <- rbind(results, read_tsv(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Urate/PWAS_CI/chr", chr, ".out")))
    if(chr==6){
        results <- rbind(results, read_tsv(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Urate/PWAS_CI/chr", chr, ".out.MHC")))
    }
}

write_tsv(results, paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Results/Urate/PWAS_CI.out"))
