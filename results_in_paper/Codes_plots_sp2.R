
## Suppl Fig 2

library(latex2exp)
library(ggplot2)


eqtls = read_excel("ExtendedDataFig2.xlsx", sheet = "dat")
eqtls.2 = read_excel("ExtendedDataFig2.xlsx", sheet = "leg")


im1 <- ggplot(eqtls, aes(x = 1:49,y=V2, size=sample)) +
  geom_point(alpha=1,color = eqtls$cls)+  
  theme(plot.title = element_text(hjust = 0.5,size = 7),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black",size = 0.5),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text = element_text(size = 6)) +
  labs(title = "Overlap with eQTLs (GTEx V8)", x="Tissues",y="Proportion")+
  scale_x_continuous(breaks = NULL)+
  coord_cartesian(ylim = c(0,0.5)) + scale_fill_manual(values = as.character(eqtls$cls))


im2 <- ggplot(eqtls, aes(x = 1:49,y=V3, size=sample)) +
  geom_point(alpha=1,color = eqtls$cls)+ 
  theme(plot.title = element_text(hjust = 0.5,size = 7),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black",size = 0.5),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text = element_text(size = 6)) +
  labs(title = "Colocalization with eQTLs (GTEx V8)", x="Tissues",y="Proportion")+
  scale_x_continuous(breaks = NULL)+
  coord_cartesian(ylim = c(0,0.25)) + 
  scale_fill_manual(values = as.character(eqtls$cls))

# eqtls.2$tissues <- as.factor(eqtls.2$tissues)
# eqtls$tissues <- as.factor(eqtls$tissues)
myColors <- eqtls$cls
names(myColors) <- eqtls$tissues

im3 <- ggplot(eqtls.2, aes(x = 1:49,y = V3)) + 
  geom_point(aes(color = tissues)) +
  scale_color_manual(name = "GTEx V8 tissues",
                     values = myColors) +
  theme(
    legend.key = element_blank(),
    legend.key.size = unit(2, "mm"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    title = element_text(size = 7),
    text = element_text(size = 6)
  ) +
  guides(color=guide_legend(ncol = 1))

pm3 <- as_ggplot(get_legend(im3))

###############################################################
###############################################################
###############################################################


p <- ggarrange(ggarrange(im1, im2,
                         nrow = 2, labels = c("a", "b"),
                         heights = c(0.5,0.5)),
               pm3,
               ncol = 2, 
               labels = c(NA, NA),
               widths = c(0.7,0.3)
)

ggsave(filename=paste0("ExtendedDataFigure2.pdf"), 
       plot=p, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/", 
       width=180, height=120, units="mm", dpi=320)



