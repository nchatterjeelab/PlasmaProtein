
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(plink2R))


lookup <-  readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")

prot <- read_tsv("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_AA_hg38.pos")

seqid <- prot$WGT
seqid <- seqid[str_detect(seqid,"RDat")]
seqid <- gsub(".wgt.RDat","",seqid)

seqid <- intersect(seqid,readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt'))
prot <- prot[match(paste0(seqid,".wgt.RDat"),prot$WGT),]


pred_matrix = matrix(nrow=696,ncol=nrow(prot))

for (i in 1:length(seqid)) {
    chr <- prot$CHR[i]

    # load in genotype files by chromosome, restrict to matching SNPs and combine
    genos = read_plink(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/window2M/byseq/",seqid[i]),impute="avg")
    if(i==1){ sampleID <- genos$fam }
    sdtmp <- apply(genos$bed, MARGIN = 2, FUN=sd)
    a <- which(sdtmp==0)
    if(length(a)>0){
        genos$bed = scale(genos$bed[,-a])
        genos$bim = genos$bim[-a,]
    }else{
        genos$bed = scale(genos$bed)
    }
    rsid <- lookup$rsid[match(genos$bim$V2,lookup$SNPid)]
    genos$bim$V2 <- rsid

    # list of SNPs that overlap loaded features
    genos.keep = rep(F,nrow(genos$bim))

    load(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_AA/",seqid[i],".wgt.RDat"))

    wgt.matrix[is.na(wgt.matrix)] = 0
    # Match up the SNPs and weights
    m = match( snps[,2] , genos$bim[,2] )
    m.keep = !is.na(m)
    snps = snps[m.keep,]
    wgt.matrix = wgt.matrix[m.keep,]
    genos.keep[ m[m.keep] ] = T

    cur.genos = genos$bed[,m[m.keep]]
    cur.bim = genos$bim[m[m.keep],]
    # Flip WEIGHTS for mismatching alleles
    #qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
    #wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

    # Predict into reference
    pred_matrix[,i] = cur.genos %*% wgt.matrix[ , "enet" ]

    rm(list = c("genos","sdtmp","a","cur.genos","cur.bim","genos.keep"))

    print(i)
}
colnames(pred_matrix) <- seqid

save(pred_matrix,sampleID,file = "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/ARICmodel_prediction/prediction.RData")


###

library(readr)

load("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/ARICmodel_prediction/prediction.RData")

seqid <- colnames(pred_matrix)
R2_aask <- numeric()
p_aask <- numeric()
R2_aric <- numeric()
p_aric <- numeric()
for (i in 1:length(seqid)){
    pheno <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/peernum/invrankpheno/50/", seqid[i], ".pheno"))
    pred <- pred_matrix[match(pheno$V2, sampleID$V2),i]
    fit <- lm(pheno$V3~pred)
    reg <- summary(fit)
    R2_aask[i] <- reg$adj.r.sq
    if("pred" %in% rownames(reg$coef)){ p_aask[i] <- reg$coef["pred","Pr(>|t|)"] }else{ p_aask[i] <- 1 }

    load(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_AA/",seqid[i],".wgt.RDat"))
    R2_aric[i] <- hsq[1]
    p_aric[i] <- cv.performance["pval","enet"]

    print(i)
}
R2 <- data.frame(seqid,R2_aask,p_aask,R2_aric,p_aric,stringsAsFactors = F)
write_tsv(R2, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/ARICmodel_prediction/R2.txt")

fdr_aask <- p.adjust(R2$p_aask, method = "fdr") # 0.789272
#fdr_aric <- p.adjust(R2$p_aric, method = "fdr")


