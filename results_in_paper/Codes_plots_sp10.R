
## Suppl Fig 8

library(gaston)
library(ggplot2)
library(GGally)


pc2 <- read.table("ExtendedDataFig10.txt",heade = T)

pc.pr <- ggpairs(pc2[,3:7],aes(color = pc2$V1),upper = list(continuous = "points"),diag = list(continuous = "blank")) + 
  scale_color_manual(values=c("#238b45","#2171b5")) + 
  theme(panel.background = element_blank(),legend.position = "bottom",axis.text = element_text(size = 6),axis.title = element_text(size = 7)) + 
  scale_x_continuous(breaks = c(-0.01,0,0.01),labels = c(-0.01,0,0.01)) + 
  scale_y_continuous(breaks = c(-0.01,0,0.01),labels = c(-0.01,0,0.01))

ggsave(filename=paste0("ExtendedDataFigure10.pdf"), 
     plot=pc.pr, device="pdf",
     path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/", 
     width=180, height=180, units="mm", dpi=320)



