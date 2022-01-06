
## Fig 1

library(ggplot2)

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 7),
  text = element_text(size = 6)
  # axis.title.x = element_text(size = 10),
  # axis.text.x = element_text(size = 8),
  # axis.title.y = element_text(size = 10),
  # axis.text.y = element_text(size = 8),
  # legend.title = element_text(size = 10)
  # legend.text = element_text(size = 8)
)


###############################################################
###############################################################
###############################################################

#### A

library(readxl)

df.peer <- read_excel("Fig1.xlsx", sheet = "1a")
df.peer_highlight <- df.peer[c(9,21),]

p1 <- ggplot(data = df.peer, aes(x = NPEER, y=NpGene, col=Race)) + 
  geom_line() + geom_point() +
  geom_point(data=df.peer_highlight,aes(x = NPEER, y=NpGene, color=Race, fill=Race), size=4,shape=18) +
  My_Theme +
  labs(y = "# significant SOMAmers      ", x="# PEER factors", title= NULL) +
  scale_colour_manual(values=c("#238b45","#2171b5")) + 
  scale_fill_manual(values=c("#238b45","#2171b5"))

###############################################################
###############################################################
###############################################################

##### B
library(ggVennDiagram)

library(limma)
library(tidyverse)
library(ggforce)

df.venn <- read_excel("Fig1.xlsx", sheet = "1b")

p2 <-  ggplot(df.venn, aes(x0 = x, y0 = y, r = r, fill = labels)) +
  geom_circle(alpha = .4, size = 1, colour = NA) +
  theme_void() +
  theme(legend.position="none")+
  My_Theme+
  annotate("text", x = -0.34, y =0.4, label = "1,618 in AA", size=2)+
  annotate("text", x = 0, y = 0, label = "1,447 overlapping", size=2)+
  annotate("text", x = 0.4, y =0.43, label = "2,004 in EA", size=2)+
  coord_fixed()+
  scale_fill_manual(values=c("#238b45","#2171b5"))

###############################################################
###############################################################
###############################################################

## C

df.effmaf <- read_excel("Fig1.xlsx", sheet = "1c")

p3 <- ggplot(data = df.effmaf, aes(x = maf, y = abs(BETA), col=Race)) +
  geom_point(alpha=0.5,size=0.5) +
  # geom_point(alpha=0.5,size=power) +  
  geom_smooth(method="lm", mapping = aes(weight = 1/power), col="#fd8d3c") +
  geom_smooth(method="lm", col="#525252") +
  scale_x_continuous(breaks=c(0, 0.1, 0.2)) +
  theme(legend.position="none")+
  My_Theme+
  facet_wrap(~Race,  ncol=2)+
  labs(x="MAF(1-MAF)", y = "Effect size", title = NULL)+
  scale_colour_manual(values=c("#238b45","#2171b5"))

###############################################################
###############################################################
###############################################################

## D

df.efftss <- read_excel("Fig1.xlsx", sheet = "1d")

p4 <- ggplot(data = df.efftss, aes(x = dist, y = abs(BETA), col=Race)) + 
  geom_point(alpha=0.5,size=0.5) +  
  scale_x_continuous(breaks=c(-0.25*10^6,0.25*10^6)) +
  theme(legend.position="none")+
  My_Theme+
  facet_wrap(~Race,  ncol=2)+
  labs(x="Distance to TSS", y = "Effect size", title = NULL)+
  scale_colour_manual(values=c("#238b45","#2171b5"))+
  coord_cartesian(ylim = c(0.2,NA)) 

###############################################################
###############################################################
###############################################################

## E

df.cond <- read_excel("Fig1.xlsx", sheet = "1e")

p5 <- ggplot(df.cond, aes(x = count)) + 
  geom_histogram(binwidth = 1, aes(fill=Race), col="white",alpha=0.8) +
  theme(legend.position="none")+
  My_Theme+
  labs(x = "# conditionally significant cis-SNPs", 
       y = "# significant SOMAmers",title=NULL) + 
  facet_grid(cols = vars(Race))+
  scale_fill_manual(values=c("#238b45","#2171b5"))


###############################################################
###############################################################
###############################################################

## Put them together

library(ggpubr)
p <- ggarrange(ggarrange(p1, p2,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.6,0.4)),
               ggarrange(p3, p4, p5,
                         ncol = 3, labels = c("c", "d","e"),
                         widths = c(0.35,0.35,0.3)),
               nrow = 2, heights = c(0.4,0.6))

ggsave(filename="Figure1.pdf", 
       plot=p, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/", 
       width=180, height=105, units="mm", dpi=320)


