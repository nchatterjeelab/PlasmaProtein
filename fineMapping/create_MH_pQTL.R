rgs = commandArgs(trailingOnly=TRUE)
pname = as.character(args[1])

print("Reading European Summary file...")
sumstat.ea <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",pname,".PHENO1.glm.linear"),header = F)
n.snp <- nrow(sumstat.ea)
ss.ea <- data.frame(sumstat.ea[,c(1:3)],sumstat.ea$V9/sumstat.ea$V10)
af.ea <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/summary_stat/",pname,".afreq"),header = F)
colnames(ss.ea) <- c("CHR","POS","ID","ZSCORE.P1")

print("Reading European LD file...")
ld.ea <- readBin(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/fine-mapping/LD/",pname,".ld.bin"), what="numeric", size=4, n=(n.snp)^2)
ld.ea = matrix(ld.ea,ncol = n.snp,byrow = T)

print("Reading African Summary file...")
sumstat.aa <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/",pname,".PHENO1.glm.linear"),header = F)
n.snp <- nrow(sumstat.aa)
ss.aa <- data.frame(sumstat.aa[,c(1:3)],sumstat.aa$V9/sumstat.aa$V10)
af.aa <- read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/summary_stat/",pname,".afreq"),header = F)
colnames(ss.aa) <- c("CHR","POS","ID","ZSCORE.P1")

print("Reading African LD file...")
ld.aa <- readBin(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/fine-mapping/LD/",pname,".ld.bin"), what="numeric", size=4, n=(n.snp)^2)
ld.aa = matrix(ld.aa,ncol = n.snp,byrow = T)


indx.col = cbind(c("dodgerblue4","deepskyblue2","darkolivegreen4","goldenrod2","firebrick2"),c("[0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))

### EA:

sentinel.ea <- which.min(sumstat.ea$V12)
ld.sent.ea = ld.ea[sentinel.ea,]
range = c(sumstat.ea[sentinel.ea,2] - 500000,sumstat.ea[sentinel.ea,2] + 500000)
sub.indx = intersect(which(sumstat.ea$V2 < range[2]),which(sumstat.ea$V2 > range[1]))
ss.ea = sumstat.ea[sub.indx,]
ld.sent.ea.sub = ld.sent.ea[sub.indx]

ld.sent.ea.sub.1 = cut(ld.sent.ea.sub^2,breaks = c(0,0.2,0.4,0.6,0.8,1),include.lowest=T)
ld.sent.ea.sub.1 = unlist(lapply(ld.sent.ea.sub.1,FUN = function(x){s1 = match(x,indx.col[,2]);return(indx.col[s1,1])}))

print("reading EA fine mapping files...")
t.ea = readRDS(paste0("/dcl01/chatterj/data/diptavo/pQTL/susie/",pname,"/White/susie_10.rds"))
pip.ea = t.ea[[3]]$pip
pip.ea.sub = pip.ea[sub.indx]

print("plotting EA...")
pdf(paste0("/dcl01/chatterj/data/diptavo/pQTL/plots/mh/",pname,"_EA.pdf"))
plot(1:dim(ss.ea)[1],-log10(ss.ea$V12),pch = 19,cex = 0.8,col = ld.sent.ea.sub.1,ylab = "-log10(pvalue)",xlab = paste("Position on chromosome:",ss.ea[1,1],"(Mb)"),xaxt = "n",main = paste(pname, ": European Americans"))
axis(1,at = c(1,dim(ss.ea)[1],sentinel.ea),labels = round(c(ss.ea[1,2],ss.ea[dim(ss.ea)[1],2],ss.ea[sentinel.ea,2])/1000000,3))
dev.off()

pdf(paste0("/dcl01/chatterj/data/diptavo/pQTL/plots/fm/",pname,"_finemapping_EA.pdf"))
plot(1:dim(ss.ea)[1],pip.ea.sub,pch = 19,cex = 0.8,col = ld.sent.ea.sub.1,ylab = "Posterior Inclusion Probability",
     xlab = paste("Position on chromosome:",ss.ea[1,1],"(Mb)"),xaxt = "n",main = c(paste0("Fine Mapping: ",pname), paste0("European Americans")),ylim = c(0,1))
axis(1,at = c(1,dim(ss.ea)[1],sentinel.ea),labels = round(c(ss.ea[1,2],ss.ea[dim(ss.ea)[1],2],ss.ea[sentinel.ea,2])/1000000,3))
abline(h = 0.85,col = "Black",lwd = 2,lty = 1)
abline(h = 0.5,col = "grey",lwd = 2,lty = 2)
dev.off()

### AA

sentinel.aa <- which.min(sumstat.aa$V12)
ld.sent.aa = ld.aa[sentinel.aa,]
range = c(sumstat.aa[sentinel.aa,2] - 500000,sumstat.aa[sentinel.aa,2] + 500000)
sub.indx = intersect(which(sumstat.aa$V2 < range[2]),which(sumstat.aa$V2 > range[1]))
ss.aa = sumstat.aa[sub.indx,]
ld.sent.aa.sub = ld.sent.aa[sub.indx]

ld.sent.aa.sub.1 = cut(ld.sent.aa.sub^2,breaks = c(0,0.2,0.4,0.6,0.8,1),include.lowest=T)
ld.sent.aa.sub.1 = unlist(lapply(ld.sent.aa.sub.1,FUN = function(x){s1 = match(x,indx.col[,2]);return(indx.col[s1,1])}))

print("reading AA fine mapping files...")
t.aa = readRDS(paste0("/dcl01/chatterj/data/diptavo/pQTL/susie/",pname,"/Black/susie_10.rds"))
pip.aa = t.aa[[3]]$pip
pip.aa.sub = pip.aa[sub.indx]

print("plotting AA...")
pdf(paste0("/dcl01/chatterj/data/diptavo/pQTL/plots/mh/",pname,"_AA.pdf"))
plot(1:dim(ss.aa)[1],-log10(ss.aa$V12),pch = 19,cex = 0.8,col = ld.sent.aa.sub.1,ylab = "-log10(pvalue)",xlab = paste("Position on chromosome:",ss.aa[1,1],"(Mb)"),xaxt = "n",main = paste(pname, ": African Americans"))
axis(1,at = c(1,dim(ss.aa)[1],sentinel.aa),labels = round(c(ss.aa[1,2],ss.aa[dim(ss.aa)[1],2],ss.aa[sentinel.aa,2])/1000000,3))
dev.off()

pdf(paste0("/dcl01/chatterj/data/diptavo/pQTL/plots/fm/",pname,"_finemapping_AA.pdf"))
plot(1:dim(ss.aa)[1],pip.aa.sub,pch = 19,cex = 0.8,col = ld.sent.aa.sub.1,ylab = "Posterior Inclusion Probability",
     xlab = paste("Position on chromosome:",ss.aa[1,1],"(Mb)"),xaxt = "n",main = c(paste0("Fine Mapping: ",pname), paste0("African Americans")),ylim = c(0,1))
axis(1,at = c(1,dim(ss.aa)[1],sentinel.aa),labels = round(c(ss.aa[1,2],ss.aa[dim(ss.aa)[1],2],ss.aa[sentinel.aa,2])/1000000,3))
abline(h = 0.85,col = "Black",lwd = 2,lty = 1)
abline(h = 0.5,col = "grey",lwd = 2,lty = 2)
dev.off()
