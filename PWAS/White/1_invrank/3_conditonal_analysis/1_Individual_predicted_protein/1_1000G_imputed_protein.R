
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(plink2R))

annota <- read_tsv("/home-2/jzhan218@jhu.edu/data/jzhan218/PWAS/Results/prot.anno_autosomal.txt")
prot <- read_tsv("/home-2/jzhan218@jhu.edu/data/jzhan218/PWAS/Results/White/Plasma_Protein_0.05.pos")

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

snp[["flip"]] = a1 != ref1

return(snp)
}

pred_matrix = matrix(nrow=498,ncol=nrow(prot))

for (i in 1:nrow(prot)) {
    chr <- prot$CHR[i]

    # load in genotype files by chromosome, restrict to matching SNPs and combine
    genos = read_plink(paste0("/home-2/jzhan218@jhu.edu/data/jzhan218/PWAS/Results/LDref/EUR/chr",chr),impute="avg")
    sdtmp <- apply(genos$bed, MARGIN = 2, FUN=sd)
    a <- which(sdtmp==0)
    genos$bed = scale(genos$bed[,-a])
    genos$bim = genos$bim[-a,]

    # list of SNPs that overlap loaded features
    genos.keep = rep(F,nrow(genos$bim))

    load(paste0("/home-2/jzhan218@jhu.edu/data/jzhan218/PWAS/Results/White/coefs_remove_ambiguous_snp_0.05/",prot$WGT[i]))

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
    qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
    wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

    # Predict into reference
    pred_matrix[,i] = cur.genos %*% wgt.matrix[ , "enet" ]

    rm(list = c("genos","sdtmp","a","cur.genos","cur.bim","genos.keep"))

    print(i)
}
colnames(pred_matrix) <- unlist(lapply(prot$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)}))

write_tsv(as.data.frame(pred_matrix), paste0("/home-2/jzhan218@jhu.edu/data/jzhan218/PWAS/Results/White/1000G_imputed/1000G_imputed_prot.txt"))

