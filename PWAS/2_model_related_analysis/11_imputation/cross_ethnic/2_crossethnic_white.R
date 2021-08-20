
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
library(dplyr)
library(stringr)
library(genio)

allele.qc = function(a1,a2,ref1,ref2) {
        a1 = toupper(a1)
        a2 = toupper(a2)
        ref1 = toupper(ref1)
        ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
} # a1: minor(alt); a2: major(ref)

annota.pwas <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_AA_hg38.pos")
annota.pwas$seq <- gsub("[.].*","",annota.pwas$WGT)

RDat <- dir("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA")
RDat <- RDat[str_detect(RDat,"RDat")]
seq.w <- unlist(lapply(str_split(RDat, "[./]"), FUN=function (x){x[1]}))
RDat <- dir("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_AA")
RDat <- RDat[str_detect(RDat,"RDat")]
seq.b <- unlist(lapply(str_split(RDat, "[./]"), FUN=function (x){x[1]}))

seq <- intersect(seq.b, seq.w)

## apply EA to AA
rsq <- numeric()
rsq.cross <- numeric()
hsq.all <- numeric()
for (i in 1:length(seq)){
  gene <- seq[i]

  load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_EA/", gene, ".wgt.RDat"))

  if(sum(wgt.matrix[,"enet"]!=0)==0){
    load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_AA/", gene, ".wgt.RDat"))
    hsq.all[i] <- hsq[1]
    rsq.cross[i] <- 0
    rsq[i] <- cv.performance["rsq","enet"]
    next
  }

  wgt.matrix <- wgt.matrix[wgt.matrix[,"enet"]!=0,]

  a <- read_plink(paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/', gene))
  b <- t(a$X)

  snpid <- intersect(rownames(wgt.matrix),colnames(b))

  if(length(snpid)==0){
    rsq.cross[i] <- 0
  } else if(length(snpid)==1){
    wgt.matrix <- wgt.matrix[snpid,]
    b <- b[,snpid]
    b[is.na(b)] <- mean(b, na.rm = T)

    imputed <- b * wgt.matrix["enet"]
    tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/80/",gene,".pheno"))
    true <- tmp$V3; names(true) <- tmp$V2
    sampleid <- intersect(names(imputed), names(true))

    reg <- summary(lm( true[sampleid] ~ imputed[sampleid] ))
    rsq.cross[i] <- reg$adj.r.sq
  } else {
    wgt.matrix <- wgt.matrix[snpid,]
    b <- b[,snpid]
    b <- apply(b, MARGIN=2, FUN=function(x){x[is.na(x)] <- mean(x,na.rm=T); x})

    imputed <- b %*% matrix(wgt.matrix[,"enet"], ncol=1)
    tmp <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/80/",gene,".pheno"))
    true <- tmp$V3; names(true) <- tmp$V2
    sampleid <- intersect(rownames(imputed), names(true))

    reg <- summary(lm( true[sampleid] ~ imputed[sampleid,] ))
    rsq.cross[i] <- reg$adj.r.sq
  }

  load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_weights_AA/", gene, ".wgt.RDat"))
  hsq.all[i] <- hsq[1]
  rsq[i] <- cv.performance["rsq","enet"]
  print(i)
}
res <- data.frame(SeqId = seq, hsq = hsq.all, rsq = rsq, rsq.cross = rsq.cross)
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/CrossEthnic")
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/CrossEthnic/EAtoAA_Predicted_Accuracy.txt"))
