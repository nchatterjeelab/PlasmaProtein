
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(dplyr)
library(stringr)
library(genio)
library( snpStats )


args <- commandArgs(T)

for(i in 1:length(args)){
     eval(parse(text=args[[i]]))
}

annota.twas <- read_tsv(paste0('/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/',tissue,'.P01.pos'))

annota.pwas <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/Plasma_Protein.pos")

geneinboth <- intersect(annota.pwas$ID, annota.twas$ID)

#### correlation of imputed protein vs imputed gene expression

allele.qc <- function(a1, a2, ref1, ref2) {
  a1 = toupper(a1); a2 = toupper(a2)
  ref1 = toupper(ref1); ref2 = toupper(ref2)

	ref = ref1; flip = ref
	flip[ref == "A"] = "T"; flip[ref == "T"] = "A"; flip[ref == "G"] = "C"; flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2; flip = ref
	flip[ref == "A"] = "T"; flip[ref == "T"] = "A"; flip[ref == "G"] = "C"; flip[ref == "C"] = "G"
	flip2 = flip;

	# snp = list()
	# snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C")); snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F; snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	# snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1) | paste0(a1,a2,ref1,ref2) %in% c('ACGT','AGCT','TCGA','TGCA','CATG','CTAG','GATC','GTAC'
  snp = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}  # a1: minor(alt); a2: major(ref)

a <- read.table("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/LDREF/1000G.EUR.10.fam")
b <- read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/LDref/EUR/chr10.fam")
inds <- intersect(a$V2, b$V2)

res <- data.frame(seqID ="a", gene="a",cor=1, stringsAsFactors = F)

for (i in 1:length(geneinboth)){
  gene <- geneinboth[i]

  seq.wgt <- annota.pwas$WGT[annota.pwas$ID==gene]
  seq <- unlist(lapply(str_split(seq.wgt, "[./]"), FUN=function (x){x[1]}))
  seq <- seq[!(is.na(seq))]
  chr <- annota.pwas$CHR[annota.pwas$ID==gene][1]

  if(length(seq) == 1){

    load(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp/", seq.wgt))

    ## NOTE: read.plink = 2 - read_plink!!!
    a <- read.plink(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/LDref/EUR/chr',chr), select.snps = rownames(wgt.matrix))
    b <- as(a$genotypes, Class="numeric")
    b <- b[,rownames(wgt.matrix)]
    if(sum(is.na(b))!=0){
      b <- apply(b, MARGIN=2, FUN=function(x){x[is.na(x)] <- mean(x,na.rm=T); x})
    }

    snps1000G <- a$map[rownames(wgt.matrix),]

    qc <- allele.qc( snps$V5 , snps$V6, snps1000G$allele.1 , snps1000G$allele.2 )
    if(sum(qc) > 0){
       wgt.matrix[,"enet"][ qc ] <- -1 * wgt.matrix[,"enet"][ qc ]
    }

    impute.pwas <- b %*% matrix(wgt.matrix[,"enet"], ncol=1)
    impute.pwas <- impute.pwas[inds,]

    load(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/", annota.twas$WGT[annota.twas$ID==gene]))
    wgt.matrix <- data.frame(wgt.matrix, snps)
    wgt.matrix <- wgt.matrix[unique(rownames(wgt.matrix)),]

    a <- read.plink(paste0('/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/LDREF/1000G.EUR.',chr), select.snps = unique(wgt.matrix$V2))
    b <- as(a$genotypes, Class="numeric")
    b <- b[,wgt.matrix$V2]
    if(sum(is.na(b))!=0){
      b <- apply(b, MARGIN=2, FUN=function(x){x[is.na(x)] <- mean(x,na.rm=T); x})
    }

    snps1000G <- a$map[wgt.matrix$V2,]

    qc <- allele.qc( snps$V5 , snps$V6, snps1000G$allele.1 , snps1000G$allele.2 )
    if(sum(qc) > 0){
       wgt.matrix$enet[ qc ] <- -1 * wgt.matrix$enet[ qc ]
    }

    impute.twas <- b %*% matrix(wgt.matrix$enet, ncol=1)
    impute.twas <- impute.twas[inds,]

    tmp <- data.frame(seqID =seq, gene=gene,cor=as.numeric(cor(impute.twas,impute.pwas)), stringsAsFactors = F)
    res <- rbind(res, tmp)

  }else{
    for (j in 1:length(seq)){

      load(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp/", seq.wgt[j]))

      a <- read.plink(paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/LDref/EUR/chr',chr), select.snps = rownames(wgt.matrix))
      b <- as(a$genotypes, Class="numeric")
      b <- b[,rownames(wgt.matrix)]
      if(sum(is.na(b))!=0){
        b <- apply(b, MARGIN=2, FUN=function(x){x[is.na(x)] <- mean(x,na.rm=T); x})
      }

      snps1000G <- a$map[rownames(wgt.matrix),]

      qc <- allele.qc( snps$V5 , snps$V6, snps1000G$allele.1 , snps1000G$allele.2 )
      if(sum(qc) > 0){
         wgt.matrix[,"enet"][ qc ] <- -1 * wgt.matrix[,"enet"][ qc ]
      }

      impute.pwas <- b %*% matrix(wgt.matrix[,"enet"], ncol=1)
      impute.pwas <- impute.pwas[inds,]


      load(paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS/", annota.twas$WGT[annota.twas$ID==gene]))
      wgt.matrix <- data.frame(wgt.matrix, snps)
      wgt.matrix <- wgt.matrix[unique(rownames(wgt.matrix)),]

      a <- read.plink(paste0('/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/LDREF/1000G.EUR.',chr), select.snps = unique(wgt.matrix$V2))
      b <- as(a$genotypes, Class="numeric")
      b <- b[,wgt.matrix$V2]
      if(sum(is.na(b))!=0){
        b <- apply(b, MARGIN=2, FUN=function(x){x[is.na(x)] <- mean(x,na.rm=T); x})
      }

      snps1000G <- a$map[wgt.matrix$V2,]

      qc <- allele.qc( snps$V5 , snps$V6, snps1000G$allele.1 , snps1000G$allele.2 )
      if(sum(qc) > 0){
         wgt.matrix$enet[ qc ] <- -1 * wgt.matrix$enet[ qc ]
      }

      impute.twas <- b %*% matrix(wgt.matrix$enet, ncol=1)
      impute.twas <- impute.twas[inds,]

      tmp <- data.frame(seqID =seq.wgt[j], gene=gene,cor=as.numeric(cor(impute.twas,impute.pwas)), stringsAsFactors = F)
      res <- rbind(res, tmp)

    }
  }
  print(i)
}

res <- res[-1,]
annota.pwas <- read_tsv('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt')
res <- inner_join(res, annota.pwas[,c(1,8,9)], by=c("seqID"="seqid_in_sample"))
write_tsv(res, paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/4_CompareTWASmodels/", tissue,".txt"))

