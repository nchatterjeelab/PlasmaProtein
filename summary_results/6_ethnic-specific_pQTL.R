library(ggplot2)
## EA's pQTLs that are rare in AA
load("/Users/jnz/Document/JHU/Research/PWAS/Analysis/13_ethnic-specific_pQTL/results_EA.RData")
(sum(mafin1000GAA$MAF < 0.01)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1242455
(sum(mafin1000GAA$MAF < 0.005)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1222334
(sum(mafin1000GAA$MAF < 0.002)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09507042
(sum(mafin1000GAA$MAF == 0)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.03822938
## two counts
(sum(mafin1000GAA$MAF <= 2/659/2)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09507042

marginEA$A1 == marginAA$A1
tmp <- marginAA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginEA$BETA) ) # 0.9066282
cor(tmp,marginEA$BETA) # 0.9271907

df.EA <- data.frame(Beta_AA=tmp,Beta_EA=marginEA$BETA,ID=marginEA$ID,stringsAsFactors = F)

p.EA <- ggplot(data = df.EA, aes(x = Beta_EA, y = Beta_AA)) + 
  geom_point(size=0.5, col="#238b45") +
  geom_abline(intercept = 0, slope = 1, col="red") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black", size = 0.2),
        plot.title = element_text(size = 12, face = "bold")) +
  mdthemes::md_theme_classic() +
  labs(x = "Effect size (EA)", 
       y = "Effect size (AA)",
       title="**Common sentinel SNP of *cis*-pQTLs in EA**")
       # title=expression(paste("Common sentinel SNP for ", italic("cis"), "-pQTLs in EA")))


## AA's pQTLs that are rare in EA
load("/Users/jnz/Document/JHU/Research/PWAS/Analysis/13_ethnic-specific_pQTL/results_AA.RData")
(sum(mafin1000GEA$MAF < 0.01)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3641509
(sum(mafin1000GEA$MAF < 0.005)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3427673
(sum(mafin1000GEA$MAF < 0.002)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3
(sum(mafin1000GEA$MAF == 0)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.236478
## two counts
(sum(mafin1000GEA$MAF <= 2/498/2)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3257862

tmp <- marginEA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginAA$BETA) ) # 0.9495885
cor(tmp,marginAA$BETA) # 0.9521188

df.AA <- data.frame(Beta_EA=tmp,Beta_AA=marginAA$BETA,ID=marginAA$ID,stringsAsFactors = F)

p.AA <- ggplot(data = df.AA, aes(x = Beta_AA, y = Beta_EA)) + 
  geom_point(size=0.5,col="#2171b5") +
  geom_abline(intercept = 0, slope = 1, col="red") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black", size = 0.2)) +
  mdthemes::md_theme_classic() +
  labs(x = "Effect size (AA)", 
       y = "Effect size (EA)",
       title="**Common sentinel SNP of *cis*-pQTLs in AA**")
# title=expression(paste("Common sentinel SNP for ", italic("cis"), "-pQTLs in AA")))

p <- cowplot::plot_grid(p.EA, p.AA, ncol=2)
p

ggsave(filename="sp1.png",
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/",
       width=10, height=5, units="in", dpi=500)




