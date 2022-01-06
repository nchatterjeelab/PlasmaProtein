## Fig 2

library(ggplot2)
library(readxl)


###############################################################
###############################################################
###############################################################

as3 <- read_excel("Fig2.xlsx", sheet = "2a")
as3 <- as.data.frame(as3)
bx1 <- ggplot(as3, aes(x=as3[,2], y=as3[,1],fill = as3[,2])) +
  geom_boxplot(alpha = 0.8,notch = T,notchwidth = 0.5,outlier.shape = NA)+     scale_fill_manual(values=c("#238b45","#2171b5")) +
  coord_cartesian(ylim = c(0,50)) + theme(plot.title = element_text(hjust = 0.5,size = 7),
                                          panel.background = element_blank(),
                                          axis.text = element_text(vjust = 0.5, hjust = 0.5,size = 6),
                                          legend.position="none")+
  labs(title = "Variants in credible sets", x=NULL,y=NULL) + 
  scale_y_continuous(name = " ",limits = c(0,50),breaks = c(0,25,50))


as4 = read_excel("Fig2.xlsx", sheet = "2b")
as4 <- as.data.frame(as4)
bx2 <- ggplot(as4, aes(x=as4[,2], y=as4[,1],fill = as4[,2])) +
  geom_boxplot(alpha = 0.8,notch = T,notchwidth = 0.5,outlier.shape = NA)+     scale_fill_manual(values=c("#238b45","#2171b5")) +
  coord_cartesian(ylim = c(0,10)) + theme(plot.title = element_text(hjust = 0.5,size = 7),
                                          panel.background = element_blank(),
                                          axis.text = element_text(vjust = 0.5, hjust = 0.5,size = 6),
                                          legend.position="none")+
  labs(title = "Number of components", x=NULL,y=NULL)+ 
  scale_y_continuous(name = " ",limits = c(0,10),breaks = c(0,5,10))



ss.ea <- read_excel("Fig2.xlsx", sheet = "2c")
ss.aa.2 <- read_excel("Fig2.xlsx", sheet = "2e")

ss.ea.1 <- read_excel("Fig2.xlsx", sheet = "2d")
ss.aa.1 <- read_excel("Fig2.xlsx", sheet = "2f")

ss.ea$cols <- factor(ss.ea$cols, levels = c("lightgrey","deepskyblue2","darkolivegreen4","goldenrod2","firebrick2"))

sentinel.ea <- which(ss.ea$pos == 161106)
mh1 <- ggplot(ss.ea[1:1000,],aes(x=pos,y=p)) + 
  geom_point(alpha = 1,size = 1.2, aes(color = cols)) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=6),axis.title = element_text(size = 6),legend.position = c(0.75, 0.75),legend.text = element_text(size = 6),legend.title = element_text(size = 7),
        legend.background = element_blank(),legend.key=element_blank()) +
  scale_x_continuous(name = "",breaks = c(ss.ea$pos[1],ss.ea$pos[sentinel.ea],ss.ea$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.ea$pos[sentinel.ea],ss.ea$pos[1000])/1000000,3))) +
  scale_y_continuous(name = expression("-log"[10]*"(p-value)"),limits = c(0,450),breaks = c(0,200,400)) + 
  annotate("point", x = ss.ea$pos[314], y = ss.ea$p[314], colour = "darkmagenta",size = 3.5,shape = 18) + 
  scale_color_manual(name = expression("r"^2), 
                     labels = c("0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1"),
                     values = c("lightgrey","deepskyblue2","darkolivegreen4","goldenrod2","firebrick2")) + 
  annotate(x=ss.ea$pos[100], y=400, geom = "text",label="EA",col = "black",size = 2)

mh1 <- mh1+ guides(color = guide_legend(ncol = 2))

sentinel.aa <- which(ss.aa.2$pos == 161106)
mh2 <- ggplot(ss.aa.2[1:1000,],aes(x=pos,y=-log10(p))) + geom_point(alpha = 1,color = ss.aa.2$cols[1:1000],size = 1.2) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=6),axis.title = element_text(size = 6)) +
  scale_x_continuous(name = paste("Position on chromosome:",ss.ea[1,1],"(Mb)"), breaks = c(ss.aa.2$pos[1],ss.aa.2$pos[sentinel.aa],ss.aa.2$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.aa.2$pos[sentinel.aa],ss.ea$pos[1000])/1000000,3))) +
  scale_y_continuous(name = expression("-log"[10]*"(p-value)"),limits = c(0,130),breaks = c(0,100)) + 
  annotate("point", x = ss.aa.2$pos[262], y = -log10(ss.aa.2$p[262]), colour = "darkmagenta",size = 3.5,shape = 18) + 
  annotate(x=ss.aa.1$pos[100], y=100, geom = "text",label="AA",col = "black",size = 2)


sent.ea <- which(ss.ea.1$pos == 161106)
mh3 <- ggplot(ss.ea.1[1000:1,],aes(x=pos,y=(pip))) + geom_point(alpha = 1,color = ss.ea.1$cols[1000:1],size = 1.2) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=6),axis.title.y = element_text(size = 6)) +
  scale_x_continuous(name = "",breaks = c(ss.ea.1$pos[2],ss.ea.1$pos[sent.ea],ss.ea.1$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.ea.1$pos[sent.ea],ss.ea.1$pos[1000])/1000000,3))) +
  scale_y_continuous(name = "PIP",limits = c(0,1),breaks = c(0,0.5,1)) + 
  annotate("point", x = ss.ea.1$pos[sent.ea], y = ss.ea.1$pip[sent.ea], colour = "darkmagenta",size = 2.5,shape = 18) + 
  annotate(x=ss.ea.1$pos[100], y=0.9, geom = "text",label="EA",col = "black",size = 2)


sent.aa <- which(ss.aa.1$pos == 161106)
mh4 <- ggplot(ss.aa.1[2:1000,],aes(x=pos,y=(pip))) + geom_point(alpha = 1,color = ss.aa.1$cols[2:1000],size = 1.2) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=6),axis.title = element_text(size = 6)) +
  scale_x_continuous(name = paste("Position on chromosome:",ss.ea[1,1],"(Mb)"), breaks = c(ss.aa.1$pos[2],ss.aa.1$pos[sent.aa],ss.aa.1$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.aa.1$pos[sent.aa],ss.ea.1$pos[1000])/1000000,3))) +
  scale_y_continuous(name = "PIP",limits = c(0,1),breaks = c(0,0.5,1)) + 
  annotate("point", x = ss.aa.1$pos[sent.aa], y = ss.aa.1$pip[sent.aa], colour = "darkmagenta",size = 2.5,shape = 18) + 
  annotate(x=ss.aa.1$pos[100], y=0.8, geom = "text",label="AA",col = "black",size = 2)

library(ggpubr)
p <- ggarrange(ggarrange(bx1, bx2,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.5,0.5)),
               ggarrange(mh1, mh3,
                         ncol = 2, labels = c("c", "d"),
                         widths = c(0.5,0.5)),
               ggarrange(mh2, mh4,
                         ncol = 2, labels = c("e", "f"),
                         widths = c(0.5,0.5)),
               nrow = 3, heights = c(0.4,0.3,0.3))

ggsave(filename="Figure2.pdf", 
       plot=p, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/", 
       width=180, height=180, units="mm", dpi=320)
