#######################################
# split white and black

rm(list=ls())

system("mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38")
system("mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White")
system("mkdir /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black")

library(readr)
soma <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Data/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,1:31]
write_tsv(soma_note, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/soma_note.txt')
somaid <- soma_note$SampleId

geno.w <- read.table("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr1.fam", stringsAsFactors=F)$V2
# geno.b <- read.table("/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr1.fam", stringsAsFactors=F)$V2

length(intersect(geno.w,somaid)) ## 7213
# length(intersect(geno.b,somaid)) ##

soma.w <- soma[soma$SampleId %in% geno.w, ]
soma.b <- soma[soma$SampleId %in% geno.b, ]

pc.w <- soma.w[,c(2, (dim(soma.w)[2]-9):(dim(soma.w)[2]))]
pc.b <- soma.b[,c(2, (dim(soma.b)[2]-9):(dim(soma.b)[2]))]

write_tsv(pc.w, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pc.txt')
write_tsv(pc.b, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pc.txt')

write_tsv(soma.w, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/soma_visit_3_log2_SMP.txt')
write_tsv(soma.b, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/soma_visit_3_log2_SMP.txt')

write_delim(data.frame(cbind(0, soma.w$SampleId)), col_names=F, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/samples.txt")
write_delim(data.frame(cbind(0, soma.b$SampleId)), col_names=F, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/samples.txt")

dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/eachsoma')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/eachsoma')


#######################################
# protein annotation (TSS -- hg38 position)

rm(list=ls())

library(readr)
library(biomaRt)
library(dplyr)
library(stringr)

annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
annota <- annota[,c(1,2,5,12,25:27)]
annota <- annota[annota$flag2==0,]
annota <- annota[!(is.na(annota$uniprot_id)),]

sum(str_detect(annota$uniprot_id, "[|]")) #43 -- mapped to multiple targets?

## first search by uniprot_gn_id in the ensembl database
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
allgene <- unique(annota$uniprot_id)
se.dat <- getBM(attributes = c("uniprot_gn_id","chromosome_name", "transcription_start_site"),
                filters = "uniprot_gn_id",
                values = allgene,
                mart = ensembl)
se.dat <- se.dat[se.dat$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
se.dat <- se.dat %>%
  group_by(uniprot_gn_id) %>%
  summarise(chromosome_name=unique(chromosome_name),
            transcription_start_site=min(transcription_start_site))
tmp <- allgene[!(allgene %in% se.dat$uniprot_gn_id)] # proteins which have no results by uniprot_gn_id
prot.anno <- inner_join(annota[!(annota$uniprot_id %in% tmp),], se.dat, by=c("uniprot_id"="uniprot_gn_id"))
prot.anno <- prot.anno[,-8]

## second search by hgnc_symbol (for proteins which have no data in ensembl database by uniprot_gn_id)
allgene <- unique(annota$entrezgenesymbol[annota$uniprot_id %in% tmp] )
se.dat1 <- getBM(attributes = c("hgnc_symbol","chromosome_name", "transcription_start_site"),
      filters = "hgnc_symbol",
      values = allgene,
      mart = ensembl)
se.dat1 <- se.dat1[se.dat1$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
se.dat1 <- se.dat1 %>% group_by(hgnc_symbol) %>% summarise(chromosome_name=unique(chromosome_name), transcription_start_site=min(transcription_start_site))
prot.anno1 <- inner_join(annota[(annota$uniprot_id %in% tmp),], se.dat1, by=c("entrezgenesymbol"="hgnc_symbol"))
prot.anno1 <- prot.anno1[,-8]

prot.anno <- rbind(prot.anno, prot.anno1) ## all proteins which have data in ensembl database
write_tsv(prot.anno, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno.txt')


writeLines(prot.anno$seqid_in_sample, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid.txt')

seqid <- prot.anno$seqid_in_sample
chr <- prot.anno$chromosome_name

onX <- which(chr == 'X')
writeLines(seqid[onX], '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/onX.txt')
onY <- which(chr == 'Y')
writeLines(seqid[onY], '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/onY.txt')

prot.anno <- prot.anno[-c(onX, onY),]
write_tsv(prot.anno, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
writeLines(prot.anno$seqid_in_sample, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')


# ######################################
# ## generate covariates for pQTL analysis
#
library(foreign)
library(readr)
library(dplyr)
sample <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/samples.txt', stringsAsFactors=F)
colnames(sample) <- c('FID', 'IID')
aric.cov <- read.dta('/dcl01/chatterj/data/jzhang2/pwas/pipeline/ARIC/derive37.dta')
aric.cov <- aric.cov[,c('id','v3age31','racegrp','gender','v3center')] # ?center
aric.cov <- inner_join(sample, aric.cov, by=c('IID'='id'))
aric.cov$racegrp <- NULL
aric.cov$gender <- ifelse(aric.cov$gender == 'F', 2, 1) # F2; M1
write_tsv(aric.cov, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/aric.cov')

#######################################
# create eachsoma

rm(list=ls())

library(readr)

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
seqid <-  readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')

## white

soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,2]
soma <- soma[, seqid]
tmp <- cbind(data.frame(FID=0, IID=soma_note$SampleId), soma)
write_tsv(tmp, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/soma.pheno')

for(i in 1:length(seqid)){
  
  tmp <- cbind(data.frame(FID=0, IID=soma_note$SampleId), soma[,i])
  write_tsv(tmp, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pheno/', seqid[i], '.pheno') )
  
  print(i)
}

## black

soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,2]
soma <- soma[, seqid]
tmp <- cbind(data.frame(FID=0, IID=soma_note$SampleId), soma)
write_tsv(tmp, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/soma.pheno')

for(i in 1:length(seqid)){

  tmp <- cbind(data.frame(FID=0, IID=soma_note$SampleId), soma[,i])
  write_tsv(tmp, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pheno/', seqid[i], '.pheno') )

  print(i)
}



#########

seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid.txt')
soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,2]
soma <- soma[, seqid]

n_peer=160
peer <- readRDS(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PEER/factors.', n_peer, '.rds'))
colnames(peer) <- paste('peer', 1:n_peer, sep='')
soma <- as.matrix(soma)
sd(residuals(lm(soma~peer[,1])))
sd(soma)
cor <- cor(peer)


fit <- lm(soma~peer[,1])
af <- anova(fit)
afss <- af$"Sum Sq"

