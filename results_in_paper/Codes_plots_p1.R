library(ggplot2)

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 8),
  text = element_text(size = 7)
  # axis.title.x = element_text(size = 10),
  # axis.text.x = element_text(size = 8),
  # axis.title.y = element_text(size = 10),
  # axis.text.y = element_text(size = 8),
  # legend.title = element_text(size = 10)
  # legend.text = element_text(size = 8)
)

## Fig 1

#### A

library(readxl)

tmp <- read_excel("/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Suppl_tables.xlsx", sheet = "ST2 -- PEER factors", skip=3)
tmp <- tmp[1:20,]
df.peer <- data.frame(NPEER=rep(tmp[[1]],2), NpGene=c(tmp[[2]], tmp[[3]]), Race=c(rep("EA", nrow(tmp)), rep("AA", nrow(tmp))))
df.peer$NPEER <- as.integer(as.character(df.peer$NPEER))
df.peer <- df.peer[df.peer$NPEER<140,]
df.peer_highlight <- df.peer[as.character(c(10,30)),]

p1 <- ggplot(data = df.peer, aes(x = NPEER, y=NpGene, col=Race)) + 
  geom_line() + geom_point() +
  geom_point(data=df.peer_highlight,aes(x = NPEER, y=NpGene, color=Race, fill=Race), size=5,shape=18) +
  My_Theme +
  labs(y = "# significant SOMAmers      ", x="No. of PEER factors", title= NULL) +
  scale_colour_manual(values=c("#238b45","#2171b5")) + 
  scale_fill_manual(values=c("#238b45","#2171b5"))

##### B
library(ggVennDiagram)

library(limma)
library(tidyverse)
library(ggforce)

df.venn <- data.frame(x = c(-0.1,0.1),
                      y = c(0,0),
                      r = c(0.4,0.45),
                      labels = c('A', 'B'))
p2 <-  ggplot(df.venn, aes(x0 = x, y0 = y, r = r, fill = labels)) +
  geom_circle(alpha = .4, size = 1, colour = NA) +
  theme_void() +
  theme(legend.position="none")+
  My_Theme+
  annotate("text", x = -0.32, y =0.4, label = "1,618 in AA", size=3)+
  annotate("text", x = 0, y = 0, label = "1,447 overlapping", size=3)+
  annotate("text", x = 0.4, y =0.43, label = "2,004 in EA", size=3)+
  coord_fixed()+
  scale_fill_manual(values=c("#238b45","#2171b5"))

## C

sig.w <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_pQTL_summary_cleaned_4.0_EA.txt")
sig.b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_pQTL_summary_cleaned_4.0_AA.txt")
library(dplyr)
tmp <- read.table("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_permutations_all.thresholds_EA.txt")
sig.w <- inner_join(sig.w, tmp, by=c("SOMAmer"="V1"))
sig.w$freqvar <- sig.w$A1_AF*(1-sig.w$A1_AF)
sig.w$fR2 <- 2 * sig.w$freqvar * sig.w$Beta^2
sig.w$power <- pf(qf(sig.w$V2, df1 = 1, df2 = 7213-1-1, lower.tail = FALSE),
                  df1 = 1, df2 = 7213-1-1, lower.tail = FALSE,
                  ncp = 7213*sig.w$fR2/(1-sig.w$fR2))
tmp <- read.table("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/1_permutations_all.thresholds_AA.txt")
sig.b <- inner_join(sig.b, tmp, by=c("SOMAmer"="V1"))
sig.b$freqvar <- sig.b$A1_AF*(1-sig.b$A1_AF)
sig.b$fR2 <- 2 * sig.b$freqvar * sig.b$Beta^2
sig.b$power <- pf(qf(sig.b$V2, df1 = 1, df2 = 1871-1-1, lower.tail = FALSE),
                  df1 = 1, df2 = 1871-1-1, lower.tail = FALSE,
                  ncp = 1871*sig.b$fR2/(1-sig.b$fR2))

df.effmaf <- data.frame(BETA = c(sig.w$Beta, sig.b$Beta), 
                        maf = c(sig.w$freqvar, sig.b$freqvar),
                        power = c(sig.w$power, sig.b$power),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))

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


## D

df.efftss <- data.frame(BETA = c(sig.w$Beta, sig.b$Beta), 
                        dist = c(sig.w$TopSNP_Pos38-sig.w$TSS, sig.b$TopSNP_Pos38-sig.b$TSS),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))
p4 <- ggplot(data = df.efftss, aes(x = dist, y = abs(BETA), col=Race)) + 
  geom_point(alpha=0.5,size=0.5) +  
  scale_x_continuous(breaks=c(-0.25*10^6,0.25*10^6)) +
  theme(legend.position="none")+
  My_Theme+
  facet_wrap(~Race,  ncol=2)+
  labs(x="Distance to TSS", y = "Effect size", title = NULL)+
  scale_colour_manual(values=c("#238b45","#2171b5"))+
  coord_cartesian(ylim = c(0.2,NA)) 


## E

sig.w <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/2_conditional_analysis_1.0_EA.txt")
sig.b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/2_conditional_analysis_1.0_AA.txt")


df.cond <- data.frame(count=c(data.frame(sort(table(sig.w$SOMAmer),decreasing=T))$Freq, 
                              data.frame(sort(table(sig.b$SOMAmer),decreasing=T))$Freq),
                      Race=c(rep("EA", length(data.frame(sort(table(sig.w$SOMAmer),decreasing=T))$Freq)),
                             rep("AA", length(data.frame(sort(table(sig.b$SOMAmer),decreasing=T))$Freq))))
p5 <- ggplot(df.cond, aes(x = count)) + 
  geom_histogram(binwidth = 1, aes(fill=Race), col="white",alpha=0.8) +
  theme(legend.position="none")+
  My_Theme+
  labs(x = "# conditionally significant cis-SNPs", 
       y = "# significant SOMAmers",title=NULL) + 
  facet_grid(cols = vars(Race))+
  scale_fill_manual(values=c("#238b45","#2171b5"))



## Put them together

library(ggpubr)
p <- ggarrange(ggarrange(p1, p2,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.6,0.4)),
               ggarrange(p3, p4, p5,
                         ncol = 3, labels = c("c", "d","e"),
                         widths = c(0.35,0.35,0.3)),
               nrow = 2, heights = c(0.4,0.6))

ggsave(filename="p1.pdf", 
       plot=p, device="pdf",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/", 
       width=200, height=115, units="mm", dpi=320)


