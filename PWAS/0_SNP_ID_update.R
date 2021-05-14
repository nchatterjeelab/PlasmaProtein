
## wrongly matched rsid

#library(bigreadr)
#
#nowupdate <- fread2("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/SNPconvert/ARIC_to_RSID.txt",header=F)
#official <-  readRDS("/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/geno/all_snp_aaskandaric_lookup.rds")
#
#library(dplyr)
#dat <- left_join(nowupdate,official,by=c("V1"="SNPid"))
#
#head(dat)
#m <- dat$V2!=dat$rsid
#update <- dat[m,]
#library(readr)
#
#update <- update[,2:3]
#colnames(update) <- c("wrongname","rsid")
#write_tsv(update, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/SNPconvert/wrongupdate.txt")


## update correct rsid

library(readr)
update <- read_tsv("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/SNPconvert/wrongupdate.txt")

library(stringr)
dir.create("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_EA/")
a <- dir("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_EA")
a <- a[str_detect(a,"RDat")]
a <- gsub(".wgt.RDat","",a)
for (i in 1:length(a)){
    load(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_EA/",a[i],".wgt.RDat"))
    tmp <- snps$V2[snps$V2 %in% update$wrongname]
    if(length(tmp)>0){
        print(i)
        right <- update$rsid[match(tmp, update$wrongname)]
        m <- match(tmp,snps$V2)
        snps$V2[m] <- right
        rownames(wgt.matrix)[m] <- right
        save(snps,wgt.matrix,cv.performance,hsq,hsq.pv,N.tot,
             file=paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_EA/",a[i],".wgt.RDat"))
    }else{
        system(paste0("cp /dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_EA/",a[i],".wgt.RDat",
                      " /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_EA/",a[i],".wgt.RDat"))
    }
    system(paste0("cp /dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_EA/",a[i],".hsq",
                  " /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_EA/",a[i],".hsq"))
}




library(stringr)
dir.create("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_AA/")
a <- dir("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_AA")
a <- a[str_detect(a,"RDat")]
a <- gsub(".wgt.RDat","",a)
for (i in 1:length(a)){
    load(paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_AA/",a[i],".wgt.RDat"))
    tmp <- snps$V2[snps$V2 %in% update$wrongname]
    if(length(tmp)>0){
        print(i)
        right <- update$rsid[match(tmp, update$wrongname)]
        m <- match(tmp,snps$V2)
        snps$V2[m] <- right
        rownames(wgt.matrix)[m] <- right
        save(snps,wgt.matrix,cv.performance,hsq,hsq.pv,N.tot,
             file=paste0("/dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_AA/",a[i],".wgt.RDat"))
    }else{
        system(paste0("cp /dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_AA/",a[i],".wgt.RDat",
                      " /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_AA/",a[i],".wgt.RDat"))
    }
    system(paste0("cp /dcl01/chatterj/data/jzhang2/PWAS_tutorial/backup/Plasma_Protein_weights_AA/",a[i],".hsq",
                  " /dcl01/chatterj/data/jzhang2/PWAS_tutorial/Plasma_Protein_weights_AA/",a[i],".hsq"))
}

