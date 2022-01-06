
## Suppl Fig 8

library(readr)
library(gaston)

dat_disease <- read_tsv(paste0("ExtendedDataFig8.txt"))
  
dat1 <- dat_disease[dat_disease$disease == "Urate",]
dat2 <- dat_disease[dat_disease$disease == "Gout",]

pdf("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/ExtendedDataFigure8.pdf", width = 7, height = 3.5)
par(mfrow=c(1,2))
qqplot.pvalues(dat1$PWAS.P, col.CB ="#E41A1C32", col.abline = "red", col= alpha("#4292c6", 0.6), 
               pch=20, bty="n",
               font.main=1,
               cex.axis = .45,
               cex.lab = .62,
               cex.main = 0.64,
               main=paste0("QQ plot of Urate PWAS p-values"))

qqplot.pvalues(dat2$PWAS.P, col.CB ="#E41A1C32", col.abline = "red", col= alpha("#4292c6", 0.6), 
               pch=20, bty="n",
               font.main=1,
               cex.axis = .45,
               cex.lab = .62,
               cex.main = 0.64,
               main=paste0("QQ plot of Gout PWAS p-values"))
dev.off()


