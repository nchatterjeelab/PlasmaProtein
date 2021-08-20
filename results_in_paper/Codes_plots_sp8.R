library(readr)
library(gaston)

disease <- "Gout"
dat <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS_",disease,"_CI.txt"))
dat <- dat[!is.na(dat$PWAS.Z),]


ggsave(filename=paste0("sp9_",disease,".png"), 
       plot=qqplot.pvalues(dat$PWAS.P, col.CB ="#E41A1C32", col.abline = "red", col= alpha("#4292c6", 0.6), pch=20, bty="n",
                           main=paste0("QQ plot of PWAS p-values of ",disease)), 
       device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/", 
       width=5, height=5, units="in", dpi=500)

