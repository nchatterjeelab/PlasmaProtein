
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
  scale_colour_manual(values=c("#2171b5","#238b45")) +
  coord_cartesian(ylim = c(0.2,NA)) 


## Put them together

library(ggpubr)
p <- ggarrange(p4, p5, ncol = 2, labels = c("a", "b"))

ggsave(filename="sp2.png", 
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/", 
       width=10, height=5, units="in", dpi=500)


