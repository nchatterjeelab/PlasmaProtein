#######################################
# split white and black

rm(list=ls())

library(readr)
soma <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Data/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,1:31]
write_tsv(soma_note, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/soma_note.txt')
somaid <- soma_note$SampleId

geno.w <- read.table('/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/white/chr1.fam',stringsAsFactors=F)$V1
remove.w <- read_csv('/dcs01/arking/ARIC_static/ARIC_Data/GWAS/White/6_Documentation/f3v2_indiv_1106.csv')
remove.w <- remove.w$IID[!is.na(remove.w$recommend_drop)]
geno.w <- geno.w[!(geno.w %in% remove.w)]


geno.b <- read.table('/dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/black/chr1.fam',stringsAsFactors=F)$V1
remove.b <- read_csv('/dcs01/arking/ARIC_static/ARIC_Data/GWAS/Black/6_Documentation/ARIC_B_f3v2_indiv.csv')
remove.b <- remove.b$IID[remove.b$RecommendDrop2==1]
geno.b <- geno.b[!(geno.b %in% remove.b)]

length(intersect(geno.w,somaid)) ## 7339
length(intersect(geno.b,somaid)) ## 1871

soma.w <- soma[soma$SampleId %in% geno.w, ]
soma.b <- soma[soma$SampleId %in% geno.b, ]

pc.w <- soma.w[,c(2, (dim(soma.w)[2]-9):(dim(soma.w)[2]))]
pc.b <- soma.b[,c(2, (dim(soma.b)[2]-9):(dim(soma.b)[2]))]

write_tsv(pc.w, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/pc.txt')
write_tsv(pc.b, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/pc.txt')

write_tsv(soma.w, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/soma_visit_3_log2_SMP.txt')
write_tsv(soma.b, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/soma_visit_3_log2_SMP.txt')

write_delim(data.frame(cbind(soma.w$SampleId, soma.w$SampleId)), col_names=F, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/samples.txt")
write_delim(data.frame(cbind(soma.b$SampleId, soma.b$SampleId)), col_names=F, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/samples.txt")

dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/eachsoma')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/eachsoma')


######################################
protein annotation

rm(list=ls())

library(readr)


annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
annota <- annota[,c(1,2,5,12,25:27)]
annota <- annota[annota$flag2==0,]
annota <- annota[!(is.na(annota$uniprot_id)),]

# annota <- as.matrix(annota)
# annota[annota[,"entrezgenesymbol"]=="APOE",]
# print(annota[annota$seqid_in_sample=='SeqId_14146_92',],width=100)
# print(prot.anno[prot.anno$seqid_in_sample=='SeqId_14146_92',],width=100)
# se.dat[se.dat$uniprot_gn_id=='P68431',]
# se.dat1[se.dat1$hgnc_symbol=='HIST1H3A',]



library(biomaRt)
library(dplyr)



library(stringr)
sum(str_detect(annota$uniprot_id, "[|]"))
str_split(annota$uniprot_id, "[|]")
annota$uniprot_id

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
allgene <- unique(annota$uniprot_id)
se.dat <- getBM(attributes = c("uniprot_gn_id","chromosome_name",
                     "start_position", "end_position"),
      filters = "uniprot_gn_id",
      values = allgene,
      mart = ensembl)
se.dat <- se.dat[se.dat$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
se.dat <- se.dat %>% group_by(uniprot_gn_id) %>% summarise(chromosome_name=unique(chromosome_name), start_position=min(start_position), end_position=max(end_position))
tmp <- allgene[!(allgene %in% se.dat$uniprot_gn_id)]
se.dat$uniprot_id <- se.dat$uniprot_gn_id
prot.anno <- inner_join(annota[!(annota$uniprot_id %in% tmp),], se.dat, by="uniprot_id")
prot.anno <- prot.anno[,-8]

allgene <- unique(annota$entrezgenesymbol[annota$uniprot_id %in% tmp] )
se.dat1 <- getBM(attributes = c("hgnc_symbol","chromosome_name",
                     "start_position", "end_position"),
      filters = "hgnc_symbol",
      values = allgene,
      mart = ensembl)
se.dat1 <- se.dat1[se.dat1$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
se.dat1$entrezgenesymbol <- se.dat1$hgnc_symbol
prot.anno1 <- inner_join(annota[(annota$uniprot_id %in% tmp),], se.dat1, by="entrezgenesymbol")
prot.anno1 <- prot.anno1[,-8]

# allgene[!(allgene %in% se.dat1$entrezgenesymbol)]

prot.anno <- rbind(prot.anno, prot.anno1)

write_tsv(prot.anno, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/prot.anno.txt')

######################################
## generate covariates for pQTL analysis

library(foreign)
library(readr)
library(dplyr)
sample <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/samples.txt', stringsAsFactors=F)
colnames(sample) <- c('id', 'IID')
aric.cov <- read.dta('/dcl01/chatterj/data/jzhang2/pwas/pipeline/ARIC/derive37.dta')
aric.cov <- aric.cov[,c('id','v3age31','racegrp','gender','v3center')] # ?center
aric.cov <- inner_join(sample, aric.cov, by='id')
colnames(aric.cov)[1] <- 'FID'
aric.cov$racegrp <- NULL
aric.cov$gender <- ifelse(aric.cov$gender == 'F', 2, 1)
write_tsv(aric.cov, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/ARIC/aric.cov')

aric.cov <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/ARIC/aric.cov', header=TRUE)
pc <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/pc.txt', header=TRUE)
pc <- pc[,-1]
for(i in c(21:25)*10){
  peer <- readRDS(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/PEER/factors.', i, '.rds'))
  colnames(peer) <- paste('peer', 1:i, sep='')
  covariates <- cbind(aric.cov, pc, peer)
  write_tsv(covariates, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/pQTL/covariates/covariates.', i, '.cov'))
}

i=500

#######################################
# create eachsoma

rm(list=ls())

library(readr)

annota <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/prot.anno.txt')
soma <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Data/soma_visit_3_log2_SMP.txt', n_max=1)
soma.id <- annota$seqid_in_sample
writeLines(soma.id, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/seqid.txt')

############
## white

############

soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,2]
soma <- soma[, soma.id]
tmp <- cbind(data.frame(FID=soma_note$SampleId, IID=soma_note$SampleId), soma)
write_tsv(tmp, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/soma.pheno')

for(i in 1:length(soma.id)){
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/eachsoma/', soma.id[i]))
  
  # tmp <- cbind(soma_note, soma[,i])
  # write_tsv(tmp, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/eachsoma/', soma.id[i], '/soma.txt') )
  
  tmp <- cbind(data.frame(FID=soma_note$SampleId, IID=soma_note$SampleId), soma[,i])
  write_tsv(tmp, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/eachsoma/', soma.id[i], '/soma.pheno') )
  
  print(i)
}

############
## black
############


soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,2]
soma <- soma[, soma.id]
tmp <- cbind(data.frame(FID=soma_note$SampleId, IID=soma_note$SampleId), soma)
write_tsv(tmp, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/soma.pheno')

soma <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/soma_visit_3_log2_SMP.txt')
soma_note <- soma[,2]
soma <- soma[, soma.id]
for(i in 1:length(soma.id)){
  dir.create(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/eachsoma/', soma.id[i]))
  
  # tmp <- cbind(soma_note, soma[,i])
  # write_tsv(tmp, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/eachsoma/', soma.id[i], '/soma.txt') )
  
  tmp <- data.frame(FID=soma_note, IID=soma_note, soma[,i])
  write_tsv(tmp, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/Black/eachsoma/', soma.id[i], '/soma.pheno') )
  print(i)
}

