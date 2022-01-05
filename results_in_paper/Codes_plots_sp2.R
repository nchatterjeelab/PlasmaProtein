eqtls <- read.table("~/ExtendedDataFig2.csv",header = T)
eqtls.2 <- read.csv("~/ExtendedDataFig2_leg.csv",header = T)



im1 <- ggplot(eqtls, aes(x = 1:49,y=V2, size=sample)) +
  geom_point(alpha=1,color = eqtls$cls)+  theme(plot.title = element_text(hjust = 0.5,size = 7),axis.title.x = element_text(size = 6),axis.title.y = element_text(size = 6),
                                                panel.background = element_blank(),axis.text.x = element_blank(),axis.line = element_line(color = "black",size = 0.5),legend.text = element_text(size = 6),
                                                legend.title = element_text(size = 6),axis.text = element_text(size = 6))+
  labs(title = "Overlap with eQTLs (GTEx V8)", x="Tissues",y="Proportion")+
  scale_x_continuous(breaks = NULL)+
  coord_cartesian(ylim = c(0,0.5)) + scale_fill_manual(values = as.character(eqtls$cls))


im2 <- ggplot(eqtls, aes(x = 1:49,y=V3, size=sample)) +
  geom_point(alpha=1,color = eqtls$cls)+  theme(plot.title = element_text(hjust = 0.5,size = 7),axis.title.x = element_text(size = 6),axis.title.y = element_text(size = 6),
                                                panel.background = element_blank(),axis.text.x = element_blank(),axis.line = element_line(color = "black",size = 0.5),legend.text = element_text(size = 6),
                                                legend.title = element_text(size = 6),axis.text = element_text(size = 6))+
  labs(title = "Colocalization with eQTLs (GTEx V8)", x="Tissues",y="Proportion")+
  scale_x_continuous(breaks = NULL)+
  coord_cartesian(ylim = c(0,0.25)) + scale_fill_manual(values = as.character(eqtls$cls))

im3 <- ggplot(eqtls.2,aes(x = 1:49,y = V3)) + geom_point(aes(color = as.character(eqtls.2$cls))) +
  scale_color_manual(drop = FALSE,name = "GTEx V8 tissues", values = as.character(eqtls.2$cls),labels = eqtls$tissues) + theme(legend.background = element_blank(),legend.text = element_text(size = 3),legend.key = element_blank(),legend.title = element_text(size = 7))
pm3 <- as_ggplot(get_legend(im3))


