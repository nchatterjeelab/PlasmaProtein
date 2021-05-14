

a <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/results")
RES <- matrix(nrow = length(a), ncol = 4)
for (i in 1:length(a)){
  load(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/results/",a[i]))
  RES[i,] <- c(mean(multi1_cv_perf), mean(singl1_cv_perf), mean(multi2_cv_perf), mean(singl2_cv_perf))
}

mean(RES[,1]>RES[,2])
mean(RES[,3]>RES[,4])

colnames(RES) <- c("EUR_cvR2(joint)","EUR_cvR2(single)", "AFR_cvR2(joint)","AFR_cvR2(single)")
saveRDS(RES, "/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/jointPWAS/RES.rds")