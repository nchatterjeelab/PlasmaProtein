args = commandArgs(trailingOnly=TRUE)
ts <- as.character(args[1])

library(coloc)

ea.prots <- readLines("/dcl01/chatterj/data/diptavo/pQTL/cis_proteins.ea.txt")
prot = read.table("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt",sep = "\t",header = T)
genes <- read.table(paste0("/dcl01/chatterj/data/diptavo/pQTL/gtex/GTEx_Analysis_v8_eQTL/",ts,".v8.egenes.txt"),header = T,sep = "\t")
prot.info = prot[match(ea.prots,prot$seqid_in_sample),]

pph4 <- NULL
system.time(for(i in 1:length(ea.prots)){
  print(i)

  pnames <- ea.prots[i]
  ensg <- genes[which(genes$gene_name == as.character(prot.info$entrezgenesymbol[i])),1]
  print(pnames)
  if(length(ensg) > 0){
    print("reading sumstat")
    sumstat = read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",pnames,".PHENO1.glm.linear"),header = F)

    print("reading eqtl")
    eqtl = system(paste0("grep -w ",ensg," /dcl01/chatterj/data/diptavo/pQTL/gtex/GTEx_Analysis_v8_eQTL/",ts,".allpairs.txt > /dcl01/chatterj/data/diptavo/pQTL/gtex/coloc/",ts,".tmp"))
    eqtl = read.table(paste0("/dcl01/chatterj/data/diptavo/pQTL/gtex/coloc/",ts,".tmp"))
    eqtl.pos = unlist(lapply(as.character(eqtl$V2),FUN = function(x){strsplit(x,split = "_")[[1]][2]}))
    common = intersect(as.character(sumstat$V2),as.character(eqtl.pos))
    eqtl = eqtl[match(as.character(common),as.character(eqtl.pos)),]
    sumstat.1 = sumstat[match(as.character(common),as.character(sumstat$V2)),]
    eqtl$V2 = sumstat.1$V3
    eqtl = eqtl[which(eqtl$V6 !=0),]
    mafs = read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",pnames,".afreq"),header = F)
    mafs = mafs[match(sumstat.1$V3,mafs$V2),]

    ds2 = list(sumstat.1$V3,rep(7213,dim(sumstat.1)[1]),sumstat.1[,c(12)],mafs$V5)
    names(ds2) = c("snp","N","pvalues","MAF")
    ds1 = list(eqtl$V2,rep(653,dim(eqtl)[1]),eqtl[,c(7)],eqtl$V6)
    names(ds1) = c("snp","N","pvalues","MAF")

    ds2$snp = as.character(ds2$snp);
    ds1$snp = as.character(ds1$snp)
    ds1$type = ds2$type = "quant"
    system.time({c1 = coloc.abf(ds1, ds2)})
    pph4 <- rbind(pph4,c1$summary[-1])

  }else{
    pph4 <- rbind(pph4,rep(NA,5))
  }

  write.table(pph4,paste0("/dcl01/chatterj/data/diptavo/pQTL/gtex/coloc/",ts,".coloc.txt"),col.names = T,row.names = F,quote = F)

})
