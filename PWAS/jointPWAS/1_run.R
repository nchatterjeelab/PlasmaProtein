

args <- commandArgs(T)

for(i in 1:length(args)){
     eval(parse(text=args[[i]]))
}

# note: x and y should be scaled
jointPWAS_cv <- function (X1, Y1, X2, Y2, shared1, shared2, nfold, nc, n1, n2, alpha1=0.5, alpha2=0.5){

  set.seed(length(Y1)+length(Y2))

  N1 <- nrow(X1); tmp <- sample(N1); X1 <- X1[ tmp , ]; Y1 <- Y1[ tmp ]
  N2 <- nrow(X2); tmp <- sample(N2); X2 <- X2[ tmp , ]; Y2 <- Y2[ tmp ]

  folds1 <- cut(seq(1,N1),breaks=nfold,labels=FALSE)
  folds2 <- cut(seq(1,N2),breaks=nfold,labels=FALSE)

  multi_pred1 <- matrix(NA, nrow = N1, ncol = n1*n2*nc)
  multi_pred2 <- matrix(NA, nrow = N2, ncol = n1*n2*nc)

  singl_pred1 <- matrix(NA, nrow = N1, ncol = n1)
  singl_pred2 <- matrix(NA, nrow = N2, ncol = n2)

  for ( i in 1:nfold ) {

    indx1 = which(folds1==i,arr.ind=TRUE)
    indx2 = which(folds2==i,arr.ind=TRUE)

    y1 = Y1[-indx1]; y2 = Y2[-indx2]
    x1 <- X1[ -indx1 ,]; x2 <- X2[ -indx2 ,]

    summ1 <- numeric(ncol(x1))
    for (i in 1:ncol(x1)){
      summ1[i] <- cor(x1[,i], y1)
    }
    summ2 <- numeric(ncol(x2))
    for (i in 1:ncol(x2)){
      summ2[i] <- cor(x2[,i], y2)
    }
    R1 <- cor(x1); R2 <- cor(x2)

    res_multi <- enet_multiethnic(summ1,  summ2, R1,  R2, shared1,  shared2, n1 , n2, nc, alpha1, alpha2)
    res_singl1 <-  enet_singlethnic( summ1, R1, n1, alpha1)
    res_singl2 <-  enet_singlethnic( summ2, R2, n2, alpha2)


    multi_pred1[indx1, ] <- X1[ indx1 ,] %*% res_multi$b1
    multi_pred2[indx2, ] <- X2[ indx2 ,] %*% res_multi$b2

    singl_pred1[indx1, ] <- X1[ indx1 ,] %*% res_singl1$b
    singl_pred2[indx2, ] <- X2[ indx2 ,] %*% res_singl2$b

  }

  alltuning <- n1*n2*nc
  multi1_cv_perf <- numeric(alltuning)
  multi2_cv_perf <- numeric(alltuning)
  singl1_cv_perf <- numeric(n1)
  singl2_cv_perf <- numeric(n2)
  for (i in 1:alltuning){
    reg = summary(lm( Y1 ~ multi_pred1[,i] ))
    multi1_cv_perf[i] <- reg$adj.r.sq
    reg = summary(lm( Y2 ~ multi_pred2[,i] ))
    multi2_cv_perf[i] <- reg$adj.r.sq
  }
  for (i in 1:n1){
    reg = summary(lm( Y1 ~ singl_pred1[,i] ))
    singl1_cv_perf[i] <- reg$adj.r.sq
  }
  for (i in 1:n2){
    reg = summary(lm( Y2 ~ singl_pred2[,i] ))
    singl2_cv_perf[i] <- reg$adj.r.sq
  }

  multi1_best <- which.max(multi1_cv_perf)
  multi2_best <- which.max(multi2_cv_perf)
  singl1_best <- which.max(singl1_cv_perf)
  singl2_best <- which.max(singl2_cv_perf)

  res <- list()
  res$best_multi1 <- c(res_multi$lambda1[multi1_best], res_multi$lambda2[multi1_best], res_multi$c[multi1_best])
  res$best_multi2 <- c(res_multi$lambda1[multi2_best], res_multi$lambda2[multi2_best], res_multi$c[multi2_best])
  res$best_singl1 <- res_singl1$lambda[singl1_best]
  res$best_singl2 <- res_singl2$lambda[singl2_best]

  res_multi1 <- enet_multiethnic_bl(summ1,  summ2, R1,  R2, shared1,  shared2,
                                   res$best_multi1[1], res$best_multi1[2], res$best_multi1[3],
                                   alpha1, alpha2)
  res_multi2 <- enet_multiethnic_bl(summ1,  summ2, R1,  R2, shared1,  shared2,
                                   res$best_multi2[1], res$best_multi2[2], res$best_multi2[3],
                                   alpha1, alpha2)
  res_singl1 <-  enet_singlethnic_bl( summ1, R1, res$best_singl1 , alpha1)
  res_singl2 <-  enet_singlethnic_bl( summ2, R2, res$best_singl2 , alpha2)

  res$res_multi1 <- res_multi1
  res$res_multi2 <- res_multi2
  res$res_singl1 <-  res_singl1
  res$res_singl2 <-  res_singl2

  return(res)

}

jointPWAS_predict <- function (jointPWAS_cv, X1, X2){

  res <- list()
  res$multi1_pred <- X1 %*% matrix(jointPWAS_cv$res_multi1$b1, ncol = 1)
  res$singl1_pred <- X1 %*% matrix(jointPWAS_cv$res_singl1$b, ncol = 1)
  res$multi2_pred <- X2 %*% matrix(jointPWAS_cv$res_multi2$b2, ncol = 1)
  res$singl2_pred <- X2 %*% matrix(jointPWAS_cv$res_singl2$b, ncol = 1)

  return(res)

}

ncv1=3
ncv2=5

library(plink2R)
library(Rcpp)
library(RcppArmadillo)


setwd("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/jointPWAS/")
sourceCpp('jointPWAS.cpp')

## EA
genos = read_plink(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq_remove_ambiguous_snp/", prot),impute="avg")
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", prot, ".pheno"), stringsAsFactors = F)
pheno[,3] = scale(pheno[,3])
m = match( paste(genos$fam[,1],genos$fam[,2]) , paste(pheno[,1],pheno[,2]) )
m.keep = !is.na(m)
genos$fam = genos$fam[m.keep,]
genos$bed = genos$bed[m.keep,]
m = m[m.keep]
pheno = pheno[m,]
# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
  genos$bed <- genos$bed[,nasnps == 0]
  genos$bim <- genos$bim[nasnps == 0,]
}
genos1 <- genos; pheno1 <- pheno

## AA
genos = read_plink(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/window1M/byseq_remove_ambiguous_snp/", prot),impute="avg")
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = read.table(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/70/", prot, ".pheno"), stringsAsFactors = F)
pheno[,3] = scale(pheno[,3])
m = match( paste(genos$fam[,1],genos$fam[,2]) , paste(pheno[,1],pheno[,2]) )
m.keep = !is.na(m)
genos$fam = genos$fam[m.keep,]
genos$bed = genos$bed[m.keep,]
m = m[m.keep]
pheno = pheno[m,]
# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
  genos$bed <- genos$bed[,nasnps == 0]
  genos$bim <- genos$bim[nasnps == 0,]
}
genos2 <- genos; pheno2 <- pheno
rm(list = c("genos","pheno","nasnps"))

snps1 <- genos1$bim$V2; snps2 <- genos2$bim$V2
shared1 <- integer(length = length(snps1))
for (i in 1:length(snps1)){
  if(snps1[i] %in% snps2)
    shared1[i] <- which(snps2==snps1[i])
  else
    shared1[i] <- 0
}
shared2 <- integer(length = length(snps2))
for (i in 1:length(snps2)){
  if(snps2[i] %in% snps1)
    shared2[i] <- which(snps1==snps2[i])
  else
    shared2[i] <- 0
}

set.seed(1)
## EA
cv.all = pheno1
N = nrow(cv.all)
cv.sample = sample(N)
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,N),breaks=ncv2,labels=FALSE)
cv.sample1 <- cv.sample; cv.all1 <- cv.all; folds1 <- folds
## AA
cv.all = pheno2
N = nrow(cv.all)
cv.sample = sample(N)
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,N),breaks=ncv2,labels=FALSE)
cv.sample2 <- cv.sample; cv.all2 <- cv.all; folds2 <- folds

rm(list = c("cv.all","cv.sample","N", "folds"))

cv.performance <- matrix(NA, ncol = 4, nrow = 2)

multi1_cv_perf <- numeric(length = length(ncv2))
multi2_cv_perf <- numeric(length = length(ncv2))
singl1_cv_perf <- numeric(length = length(ncv2))
singl2_cv_perf <- numeric(length = length(ncv2))
res <- list()
for ( i in 1:ncv2 ) {

  indx1 = which(folds1==i,arr.ind=TRUE)
  indx2 = which(folds2==i,arr.ind=TRUE)

  wgt = jointPWAS_cv( X1= genos1$bed[ cv.sample1[ -indx1 ],] , Y1 = cv.all1[-indx1,3] ,
                      X2= genos2$bed[ cv.sample2[ -indx2 ],] , Y2 = cv.all2[-indx2,3] ,
                      shared1=shared1, shared2=shared2,
                      nfold=ncv1, nc=3, n1=10, n2=10, alpha1=0.5, alpha2=0.5 )

  pred <- jointPWAS_predict(jointPWAS_cv=wgt, X1=genos1$bed[ cv.sample1[ indx1 ],], X2=genos2$bed[ cv.sample2[ indx2 ],])

  Y1 <- cv.all1[indx1,3]
  Y2 <- cv.all2[indx2,3]

  reg = summary(lm( Y1 ~ pred$multi1_pred ))
  multi1_cv_perf[i] <- reg$adj.r.sq
  reg = summary(lm( Y2 ~ pred$multi2_pred ))
  multi2_cv_perf[i] <- reg$adj.r.sq
  reg = summary(lm( Y1 ~ pred$singl1_pred ))
  singl1_cv_perf[i] <- reg$adj.r.sq
  reg = summary(lm( Y2 ~ pred$singl2_pred ))
  singl2_cv_perf[i] <- reg$adj.r.sq

  res[[i]] <- list(pred,Y1,Y2,wgt)
  print(paste0("Completed cv ", i))

}

save(res,
     multi1_cv_perf,multi2_cv_perf,singl1_cv_perf,singl2_cv_perf,
     file = paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/results/", prot,".RDat"))


