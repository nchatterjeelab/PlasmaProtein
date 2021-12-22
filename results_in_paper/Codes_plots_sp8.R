
## Suppl Fig 8

library(readr)
library(gaston)

disease <- "Urate"
dat <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_",disease,"_CI.txt"))
dat1 <- dat[!is.na(dat$PWAS.Z),]

disease <- "Gout"
dat <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_",disease,"_CI.txt"))
dat2 <- dat[!is.na(dat$PWAS.Z),]


pdf("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/sp8.pdf", width = 7.87, height = 3.94)
par(mfrow=c(1,2))
qqplot.pvalues(dat1$PWAS.P, col.CB ="#E41A1C32", col.abline = "red", col= alpha("#4292c6", 0.6), 
               pch=20, bty="n",
               font.main=1,
               cex.axis = .75,
               cex.lab = .9,
               cex.main = 1.1,
               main=paste0("QQ plot of Urate PWAS p-values"))

qqplot.pvalues(dat2$PWAS.P, col.CB ="#E41A1C32", col.abline = "red", col= alpha("#4292c6", 0.6), 
               pch=20, bty="n",
               font.main=1,
               cex.axis = .75,
               cex.lab = .9,
               cex.main = 1.1,
               main=paste0("QQ plot of Gout PWAS p-values"))
dev.off()


