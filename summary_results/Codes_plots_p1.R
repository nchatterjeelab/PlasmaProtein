library(ggplot2)

## Fig 1

#### A

library(readxl)

tmp <- read_excel("/Users/jnz/Document/JHU/Research/PWAS/PAPER/Tables_updated.xlsx")
tmp <- read_excel("/Users/jnz/Document/JHU/Research/PWAS/PAPER/Tables_updated.xlsx", sheet = "ST2 -- PEER factors", skip=4)
df.peer <- data.frame(NPEER=rep(tmp[[1]],2), NpGene=c(tmp[[2]], tmp[[3]]), Race=c(rep("AA", nrow(tmp)), rep("EA", nrow(tmp))))
df.peer <- df.peer[df.peer$NPEER<140,]
df.peer_highlight <- df.peer[as.character(c(8,34)),]

p1 <- ggplot(data = df.peer, aes(x = NPEER, y=NpGene, col=Race)) + 
    geom_line() + geom_point() +
    geom_point(data=df.peer_highlight,aes(x = NPEER, y=NpGene, color=Race, fill=Race), size=5,shape=18) +
    theme(panel.background = element_blank()) +
    labs(y = "# significant SOMAmers ", x="No. of PEER factors", title= NULL) +
    scale_colour_manual(values=c("#2171b5","#238b45")) + 
    scale_fill_manual(values=c("#2171b5","#238b45"))

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
    theme(legend.position="none",
          panel.background = element_blank()) +
    annotate("text", x = -0.4, y =0.4, label = "1605 in AA")+
    annotate("text", x = 0, y = 0, label = "1442 overlapping")+
    annotate("text", x = 0.4, y =0.43, label = "1992 in EA")+
    coord_fixed()+
    scale_fill_manual(values=c("#2171b5","#238b45"))
    

## C

sig.w <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/White_allsig_cleaned.txt", stringsAsFactors = F)
sig.b <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/Black_allsig_cleaned.txt", stringsAsFactors = F)
df.cond <- data.frame(count=c(data.frame(sort(table(sig.w$V1),decreasing=T))$Freq, 
                              data.frame(sort(table(sig.b$V1),decreasing=T))$Freq),
                      Race=c(rep("EA", length(data.frame(sort(table(sig.w$V1),decreasing=T))$Freq)),
                             rep("AA", length(data.frame(sort(table(sig.b$V1),decreasing=T))$Freq))))
p3 <- ggplot(df.cond, aes(x = count)) + 
    geom_histogram(binwidth = 1, aes(fill=Race), col="white",alpha=0.8) +
    theme(legend.position="none",
          panel.background = element_blank()) +
    labs(x = "# conditionally significant cis-SNPs", 
         y = "# significant SOMAmers",title=NULL) + 
    facet_grid(cols = vars(Race))+
    scale_fill_manual(values=c("#2171b5","#238b45"))


## D

library(dplyr)
tmp <- sig.w %>% group_by(V1) %>% summarise(Nsig=max(V12)+1)
sig.w <- inner_join(sig.w, tmp, by="V1")
tmp <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/White_permutations_all.thresholds.txt")
sig.w <- inner_join(sig.w, tmp, by="V1")
MAFsig.w <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/White_MAFsig.txt", stringsAsFactors = F)
sig.w <- inner_join(sig.w, MAFsig.w[,c(2,5)], by=c("V8"="V2"))
sig.w$freqvar <- sig.w$V5.y*(1-sig.w$V5.y)
sig.w$R2 <- 2 * sig.w$freqvar * sig.w$V18^2
sig.w$power <- pf(qf(sig.w$V2.y, df1 = 1, df2 = 7213-sig.w$Nsig-1, lower.tail = FALSE), 
                  df1 = 1, df2 = 7213-sig.w$Nsig-1, lower.tail = FALSE, 
                  ncp = 7213*sig.w$R2/(1-sig.w$R2))
# sig.w$power <- 1 - pchisq(qchisq(1-sig.w$V2.y, 1), df = 1, 
#                           ncp = 7213*sig.w$R2/(1-sig.w$R2))


tmp <- sig.b %>% group_by(V1) %>% summarise(Nsig=max(V12)+1)
sig.b <- inner_join(sig.b, tmp, by="V1")
tmp <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/Black_permutations_all.thresholds.txt")
sig.b <- inner_join(sig.b, tmp, by="V1")
MAFsig.b <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/Black_MAFsig.txt", stringsAsFactors = F)
sig.b <- inner_join(sig.b, MAFsig.b[,c(2,5)], by=c("V8"="V2"))
sig.b$freqvar <- sig.b$V5.y*(1-sig.b$V5.y)
sig.b$R2 <- 2 * sig.b$freqvar * sig.b$V18^2
sig.b$power <- pf(qf(sig.b$V2.y, df1 = 1, df2 = 1871-sig.b$Nsig-1, lower.tail = FALSE),
                  df1 = 1, df2 = 1871-sig.b$Nsig-1, lower.tail = FALSE, 
                  ncp = 1871*sig.b$R2/(1-sig.b$R2))
# sig.b$power <- 1 - pchisq(qchisq(1-sig.b$V2.y, 1), df = 1, 
#                           ncp = 1871*sig.b$R2/(1-sig.b$R2))

df.effmaf <- data.frame(BETA = c(sig.w$V18, sig.b$V18), 
                        maf = c(sig.w$freqvar, sig.b$freqvar),
                        power = c(sig.w$power, sig.b$power),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))

p4 <- ggplot(data = df.effmaf, aes(x = maf, y = abs(BETA), col=Race)) +
    geom_point(alpha=0.5,size=0.5) +
    # geom_point(alpha=0.5,size=power) +  
    geom_smooth(method="lm", mapping = aes(weight = 1/power), col="#cb181d") +
    geom_smooth(method="lm", col="black") +
    scale_x_continuous(breaks=c(0, 0.1, 0.2)) +
    theme(legend.position="none",
          panel.background = element_blank()) +
    facet_wrap(~Race,  ncol=2)+
    labs(x="MAF(1-MAF)", y = "Effect size", title = NULL)+
    scale_colour_manual(values=c("#2171b5","#238b45"))


## E

annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/prot.anno_autosomal.txt")
tss <- integer()
for (i in 1:nrow(sig.w)) {
    tss[i] <- annota$transcription_start_site[annota$seqid_in_sample == sig.w$V1[i]]
}
sig.w$TSS <- tss
tss <- integer()
for (i in 1:nrow(sig.b)) {
    tss[i] <- annota$transcription_start_site[annota$seqid_in_sample == sig.b$V1[i]]
}
sig.b$TSS <- tss
df.efftss <- data.frame(BETA = c(sig.w$V18, sig.b$V18), 
                 dist = c(sig.w$V10-sig.w$TSS, sig.b$V10-sig.b$TSS),
                 Race = c(rep("EA", nrow(sig.w)),
                          rep("AA", nrow(sig.b))))
p5 <- ggplot(data = df.efftss, aes(x = dist, y = abs(BETA), col=Race)) + 
    geom_point(alpha=0.5,size=0.5) +  
    theme(legend.position="none") +
    scale_x_continuous(breaks=c(-0.5*10^6,0.5*10^6)) +
    theme(legend.position="none",
          panel.background = element_blank()) +
    facet_wrap(~Race,  ncol=2)+
    labs(x="Distance to TSS", y = "Effect size", title = NULL)+
    scale_colour_manual(values=c("#2171b5","#238b45"))


## Put them together

library(ggpubr)
p <- ggarrange(ggarrange(p1, p2,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.6,0.4)),
               ggarrange(p3, p4, p5,
                         ncol = 3, labels = c("c", "d","e")),
               nrow = 2, heights = c(0.4,0.6))

ggsave(filename="p1.png", 
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/", 
       width=12, height=7.5, units="in", dpi=500)

########################################################################
########################################################################
## plot for sentinel SNPs (sp2)
########################################################################
########################################################################

## a


library(dplyr)
pQTL.w <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Tables/1_pQTL_summary_cleaned_1.0_EA.txt")
sig.w <- pQTL.w[,c("SOMAmer","TSS","TopSNP","Beta","TopSNP_Pos38","A1_AF")]
tmp <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/White_permutations_all.thresholds.txt")
sig.w <- inner_join(sig.w, tmp, by=c("SOMAmer"="V1"))
sig.w$freqvar <- sig.w$A1_AF*(1-sig.w$A1_AF)
sig.w$R2 <- 2 * sig.w$freqvar * sig.w$Beta^2
sig.w$power <- pf(qf(sig.w$V2, df1 = 1, df2 = 7213-1-1, lower.tail = FALSE), 
                  df1 = 1, df2 = 7213-1-1, lower.tail = FALSE, 
                  ncp = 7213*sig.w$R2/(1-sig.w$R2))

pQTL.b <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Tables/1_pQTL_summary_cleaned_1.0_AA.txt")
sig.b <- pQTL.b[,c("SOMAmer","TSS","TopSNP","Beta","TopSNP_Pos38","A1_AF")]
tmp <- read.table("/Users/jnz/Document/JHU/Research/pQTL/Analysis.new/Black_permutations_all.thresholds.txt")
sig.b <- inner_join(sig.b, tmp, by=c("SOMAmer"="V1"))
sig.b$freqvar <- sig.b$A1_AF*(1-sig.b$A1_AF)
sig.b$R2 <- 2 * sig.b$freqvar * sig.b$Beta^2
sig.b$power <- pf(qf(sig.b$V2, df1 = 1, df2 = 1871-1-1, lower.tail = FALSE),
                  df1 = 1, df2 = 1871-1-1, lower.tail = FALSE, 
                  ncp = 1871*sig.b$R2/(1-sig.b$R2))


df.effmaf <- data.frame(BETA = c(sig.w$Beta, sig.b$Beta), 
                        maf = c(sig.w$freqvar, sig.b$freqvar),
                        power = c(sig.w$power, sig.b$power),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))

p4 <- ggplot(data = df.effmaf, aes(x = maf, y = abs(BETA), col=Race)) +
  geom_point(alpha=0.5,size=0.5) +
  # geom_point(alpha=0.5,size=power) +  
  geom_smooth(method="lm", mapping = aes(weight = 1/power), col="#cb181d") +
  geom_smooth(method="lm", col="black") +
  scale_x_continuous(breaks=c(0, 0.1, 0.2)) +
  theme(legend.position="none",
        panel.background = element_blank()) +
  facet_wrap(~Race,  ncol=2)+
  labs(x="MAF(1-MAF)", y = "Effect size", title = NULL)+
  scale_colour_manual(values=c("#2171b5","#238b45"))


## b

df.efftss <- data.frame(BETA = c(sig.w$Beta, sig.b$Beta), 
                        dist = c(sig.w$TopSNP_Pos38-sig.w$TSS, sig.b$TopSNP_Pos38-sig.b$TSS),
                        Race = c(rep("EA", nrow(sig.w)),
                                 rep("AA", nrow(sig.b))))
p5 <- ggplot(data = df.efftss, aes(x = dist, y = abs(BETA), col=Race)) + 
  geom_point(alpha=0.5,size=0.5) +  
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(-0.5*10^6,0.5*10^6)) +
  theme(legend.position="none",
        panel.background = element_blank()) +
  facet_wrap(~Race,  ncol=2)+
  labs(x="Distance to TSS", y = "Effect size", title = NULL)+
  scale_colour_manual(values=c("#2171b5","#238b45"))


## Put them together

library(ggpubr)
p <- ggarrange(p4, p5, ncol = 2, labels = c("a", "b"))

ggsave(filename="sp2.png", 
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/", 
       width=10, height=5, units="in", dpi=500)


