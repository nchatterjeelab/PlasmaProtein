
########################################################################
## conditional analysis -- validated v7 hits in v8

suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("plink2R"))

opt <- list()
opt$TWAS_v8 <- "/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8.out"
opt$tissue_list_v7="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V7_tissue_list.txt"
opt$tissue_list_v8="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_tissue_list.txt"
opt$tissue_n_gene_v8="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_n_gene.rds"
opt$imputed_P="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_Plasma_Protein.txt"
opt$imputed_T_v8="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_FUSION_v8.ALL/"
opt$out_v7="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/ConditionalAnalysis_CI_v7/"
opt$out_v8="/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/ConditionalAnalysis_CI_v8_validated/"

dir.create(paste0(opt$out_v8))
dir.create(paste0(opt$out_v8, "/RDat"))
dir.create(paste0(opt$out_v8, "/Table"))

tissue_list_v7 <- readLines(opt$tissue_list_v7)
tissue_list_v8 <- readLines(opt$tissue_list_v8)

pred_prot <- suppressMessages(read_tsv(opt$imputed_P)) # load imputed cis-regulated protein levels for reference individuals

cat(paste0("* Starting to validate conditional analysis in V8 tissues --\n"))

for (tissue in tissue_list_v7){
  if( !(tissue %in% tissue_list_v8) ){next}
  pred_ge <- suppressMessages(read_tsv(paste0(opt$imputed_T_v8,"/", tissue,".txt"))) # load imputed cis-regulated gene expression levels for reference individuals

  load(paste0(opt$out_v7,"/RDat/",tissue,".RDat"))

  twas.hit_v7 <- twas.hit
  dat.twas <- suppressMessages(read_tsv(paste0(opt$TWAS_v8,"/",tissue,".out"))) # load V8 TWAS result table
  dat.twas <- dat.twas[match(twas.hit_v7,dat.twas$ID),]

  ## perform conditional anlaysis for each sentinel PWAS gene and its nearby TWAS genes
  PcT.z <- numeric(); PcT.p <- numeric()
  TcP.z <- numeric(); TcP.p <- numeric()
  twas.p <- numeric()
  twas.BETA <- numeric(); twas.SE <- numeric(); twas.CI <- character()
  dist <- integer(); corr <- numeric()
  for (i in 1:nrow(dat.sentinel.pwas)) {
    chr <- dat.sentinel.pwas$CHR[i]
    dist[i] <- dat.twas$P0[i] - dat.sentinel.pwas$P0[i]

    twas.p[i] <- dat.twas$TWAS.P[i]
    twas.hit[i] <- dat.twas$ID[i]

    twas.BETA[i] <- dat.twas$TWAS.BETA[i]
    twas.SE[i] <- dat.twas$TWAS.SE[i]
    twas.CI[i] <- dat.twas$TWAS.CI[i]

    pred_mat = matrix(nrow=498,ncol=2)
    tmp <- strsplit(dat.sentinel.pwas$FILE[i], "/")[[1]]; tmp <- tmp[length(tmp)]; tmp <- substr(tmp, start=1, stop=nchar(tmp)-9)
    pred_mat[,1] = pred_prot[[tmp]]
    tmp <- strsplit(dat.twas$FILE[i], "/")[[1]]; tmp <- tmp[length(tmp)]; tmp <- substr(tmp, start=nchar(dat.twas$PANEL[i])+2, stop=nchar(tmp)-9)
    if(is.na(tmp)) { pred_mat[,2] <- 0; corr[i] <- NA; PcT.z[i] = NA; PcT.p[i] = NA; TcP.z[i] = NA; TcP.p[i] = NA; next }
    pred_mat[,2] = pred_ge[[tmp]]

    pred_mat = scale( pred_mat )

    # compute cis-regulated genetic correlation
    corr_mat = cor(pred_mat)
    corr[i] = corr_mat[1,2]

    b.se = sqrt(1 - corr_mat[1,2]^2)

    # estimate conditional effect size (P conditional on T, P|T)
    b = dat.sentinel.pwas$PWAS.Z[i] - corr[i] * dat.twas$TWAS.Z[i]
    PcT.z[i] = b / b.se
    PcT.p[i] = 2*(pnorm( abs( PcT.z[i] ) , lower.tail=F))

    # estimate conditional effect size (T conditional on P, T|P)
    b = dat.twas$TWAS.Z[i] - corr[i] * dat.sentinel.pwas$PWAS.Z[i]
    TcP.z[i] = b / b.se
    TcP.p[i] = 2*(pnorm( abs( TcP.z[i] ) , lower.tail=F))

  }
    save(dat.sentinel.pwas,
         PcT.z, PcT.p,
         TcP.z, TcP.p,
         twas.p,
         twas.hit_v7, twas.hit,
         twas.BETA, twas.SE, twas.CI,
         dist, corr,
         file = paste0(opt$out_v8, "/RDat/",tissue,".RDat"))

    cat(paste0(tissue, " is completed. ",sum(!is.na(twas.hit))," genes are validated. \n"))

}

############################################################
## Clean results to table
############################################################

##############################
# 1. Analyze certain tissue (usually relevant tissues)

cat(paste0("* Cleaning conditional analysis results to ./Table/\n"))

for (tissue in tissue_list_v7){
  if( !(tissue %in% tissue_list_v8) ){next}
    ## load and reformat conditional analysis outputs
    load(paste0(opt$out_v8, "/RDat/",tissue,".RDat"))
    PWAS_hit=dat.sentinel.pwas$ID
    res <- list()
    for (i in 1:length(PWAS_hit)){
        TWAS_hit_v7 <- twas.hit_v7[i]
        TWAS_hit <- twas.hit[i]
        TWAS_p <- twas.p[i]
        TWAS_BETA <- twas.BETA[i]
        TWAS_SE <- twas.SE[i]
        TWAS_CI <- twas.CI[i]
        Dist_of_hits <- dist[i]
        Corr_of_hits <- corr[i]
        PcT_p <- PcT.p[i]
        TcP_p <- TcP.p[i]

        res[[PWAS_hit[i]]] <- list(TWAS_hit_v7=TWAS_hit_v7,
                                   TWAS_hit=TWAS_hit, TWAS_p=TWAS_p,
                                   TWAS_BETA = TWAS_BETA,TWAS_SE = TWAS_SE,TWAS_CI = TWAS_CI,
                                   Dist_of_hits=Dist_of_hits, Corr_of_hits=Corr_of_hits,
                                   PcT_p=PcT_p, TcP_p=TcP_p)
    }

    ## TWAS significance level corrected by number of transcripts
    p.twas <- 0.05/(readRDS(opt$tissue_n_gene_v8)[tissue])

    ## reformat to output tables (if a TWAS gene reaches significance level, will be marked with a star)
    TWAS_p <- numeric()
    TWAS_hit <- character()
    TWAS_BETA <- numeric()
    TWAS_SE <- numeric()
    TWAS_CI <- character()
    Dist_of_hits <- integer()
    Corr_of_hits <- numeric()
    PcT_p <- numeric()
    TcP_p <- numeric()
    sig <- logical()
    for (i in 1:length(res)){
        m <- (res[[i]]$TWAS_p < p.twas)
        m[is.na(m)] <- FALSE
        sig[i] <- m
        TWAS_p[i] <- res[[i]]$TWAS_p
        TWAS_hit[i] <- res[[i]]$TWAS_hit
        TWAS_BETA[i] <- res[[i]]$TWAS_BETA
        TWAS_SE[i] <- res[[i]]$TWAS_SE
        TWAS_CI[i] <- res[[i]]$TWAS_CI
        Dist_of_hits[i] <- res[[i]]$Dist_of_hits
        Corr_of_hits[i] <- res[[i]]$Corr_of_hits
        PcT_p[i] <- res[[i]]$PcT_p
        TcP_p[i] <- res[[i]]$TcP_p
    }

  cond_v7 <- read_tsv(paste0(opt$out_v7, "/Table/",tissue,".txt"), col_types = cols())

    a <- data.frame(TWAS_p=ifelse(sig,paste0(signif(TWAS_p, 3),"*"),TWAS_p),
                    TWAS_BETA=TWAS_BETA,
                    TWAS_SE=TWAS_SE,
                    TWAS_CI=TWAS_CI,

                    Corr_of_hits=signif(Corr_of_hits,3),
                    PcT_p=signif(PcT_p,3),
                    TcP_p=signif(TcP_p,3))

  colnames(a) <- paste0("v8_",colnames(a))
  a <- cbind(cond_v7, a)

    write_tsv(a, paste0(opt$out_v8, "/Table/",tissue,".txt"))
}
cat(paste0("** Completed. \n"))

##############################
# 2. All-tissue analysis

# get all the regional sentinel PWAS genes
PWAS_hit=dat.sentinel.pwas$ID
cond_v7 <- read_tsv(paste0(opt$out_v7, "/Table/all-tissue.txt"))

b <- tibble()
for (i in 1:nrow(cond_v7)){
  tissue <- gsub("[*]","",cond_v7$min_TWAS_Tissue[i])
  gene <- gsub("[*]","",cond_v7$min_TWAS_TWAS_hit[i])
  cond_v8 <- suppressMessages(read_tsv(paste0(opt$out_v8, "/Table/",tissue,".txt"), col_types = cols()))
  b <- rbind(b, cond_v8[i,15:21])
}

a <- cbind(cond_v7[,1:15],b,cond_v7[,16:17])

write_tsv(a, paste0(opt$out_v8, "/Table/all-tissue.txt"))

cat(paste0("* Validated all-tissue V8 results in ./Table/all-tissue.txt \n"))

