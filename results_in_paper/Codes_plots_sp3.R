
## Suppl Fig 3

library(latex2exp)
library(ggplot2)

ts.shr2 = read.table("ExtendedDataFig3.txt",header = F)
ts <- ggplot(ts.shr2, aes(x = ts.shr2[,1], y = ts.shr2[,2],fill = ts.shr2[,3])) + 
  geom_ribbon(stat = "smooth",aes(ymin = 0, ymax = ..y..), alpha = .5,
              method = "gam", se=FALSE, formula = y ~ s(x, k = 20)) + 
  theme(panel.background = element_blank(),axis.text = element_text(size = 7),legend.position = c(0.55, 0.45),axis.title = element_text(size = 6),legend.text = element_text(size = 6)) + scale_fill_manual(name = "", labels = c("cis-eQTLs in GTEx V8 that are also cis-pQTL","cis-eQTL in GTEx V8"), values = c("#00BFC4","#F8766D")) + xlab("Proportion") + ylab("Number of Tissues")

ggsave(filename="ExtendedDataFigure3.pdf",
       plot=ts, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/",
       width=180, height=105, units="mm", dpi=320)


