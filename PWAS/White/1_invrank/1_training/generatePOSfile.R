
rm(list=ls())

library(readr)
library(stringr)
library(dplyr)

## start and end ###
library(readr)
library(biomaRt)
annota <- read_tsv('/dcs01/arking/ARIC_static/ARIC_Data/Proteomics/ARIC-SomaLogic_Nov2019/Abbreviated annotation visits 3 and 5.txt')
annota <- annota[,c(1,2,5,12,25:27)]
annota <- annota[annota$flag2==0,]
annota <- annota[!(is.na(annota$uniprot_id)),]

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

prot.anno <- rbind(prot.anno, prot.anno1)
prot.anno <- prot.anno[,c("seqid_in_sample", "entrezgenesymbol", "chromosome_name","start_position","end_position")]

RDat <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_0.05/")
RDat <- RDat[str_detect(RDat,"RDat")]
seq <- unlist(lapply(str_split(RDat, "[.]wgt[.]RDat"), FUN=function (x){x[1]}))
prot.anno <- prot.anno[prot.anno$seqid_in_sample %in% seq,]

POS <- data.frame(PANEL="Plasma_Protein",
                  WGT=paste0(prot.anno$seqid_in_sample, ".wgt.RDat"),
                  ID=prot.anno$entrezgenesymbol,
                  CHR=prot.anno$chromosome_name,
                  P0=prot.anno$start_position,
                  P1=prot.anno$end_position,
                  N=7213)
write_tsv(POS, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/Plasma_Protein_0.05.pos")



