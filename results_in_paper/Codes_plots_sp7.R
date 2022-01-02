
## Suppl Fig 7

library(readr)
library(gaston)

dat <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_null_cont_normal.txt"))
dat <- dat[!is.na(dat$PWAS.Z),]

pdf("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/sp7.pdf", 
    width = 3.5, height = 3.5)
qqplot.pvalues(dat$PWAS.P, col.CB ="#E41A1C32", col.abline = "red", col= alpha("#4292c6", 0.4), 
               pch=20, bty="n",
               font.main=1,
               cex.axis = .45,
               cex.lab = .62,
               cex.main = 0.64,
               main=paste0("QQ plot of null PWAS p-values")
)
dev.off()

