
## Suppl Fig 1

library(latex2exp)
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

# EA (panel a)

## EA's pQTLs that are rare in AA
load("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ethnic_specific/results_EA.RData")
(sum(mafin1000GAA$MAF < 0.01)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1334333
(sum(mafin1000GAA$MAF < 0.005)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1304348
(sum(mafin1000GAA$MAF < 0.002)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09995002
(sum(mafin1000GAA$MAF == 0)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.04197901
## two counts
(sum(mafin1000GAA$MAF <= 2/659/2)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09995002

marginEA$A1 == marginAA$A1
tmp <- marginAA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginEA$BETA) ) # 0.9032445
cor(tmp,marginEA$BETA) # 0.9275709

df.EA <- data.frame(Beta_AA=tmp,Beta_EA=marginEA$BETA,ID=marginEA$ID,stringsAsFactors = F)

df.EA[(df.EA$Beta_AA > 0.2) & (df.EA$Beta_EA < -0.6),] # 865
df.EA[(df.EA$Beta_AA > 0) & (df.EA$Beta_EA < -1),] # 1391
df.EA[(df.EA$Beta_AA < 0) & (df.EA$Beta_EA > 0.5),] # 1342
df.EA[(df.EA$Beta_AA > 0.2) & (df.EA$Beta_AA < 0.3) & (df.EA$Beta_EA > 1),] # 277


df.EA_highlight <- df.EA[c(865,1391,1342,277), ]
  
p.EA <- ggplot(data = df.EA, aes(x = Beta_EA, y = Beta_AA)) + 
  geom_point(size=0.5, col="#2171b5") +
  geom_abline(intercept = 0, slope = 1, col="red") +
  theme(axis.line = element_line(color="black", size = 0.2)) +
  ylim(-2,2)+xlim(-2,2)+
  mdthemes::md_theme_classic() +
  labs(x = "Effect size (EA)", 
       y = "Effect size (AA)",
       title="Common sentinel SNP of *cis*-pQTLs in EA") +
  My_Theme +
  geom_point(data=df.EA_highlight, aes(x = Beta_EA, y = Beta_AA), size=0.5, col="darkorange")


       



###############################################################
###############################################################
###############################################################

# AA (panel b)

## AA's pQTLs that are rare in EA
load("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ethnic_specific/results_AA.RData")
(sum(mafin1000GEA$MAF < 0.01)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3665215
(sum(mafin1000GEA$MAF < 0.005)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3447418
(sum(mafin1000GEA$MAF < 0.002)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3036714
(sum(mafin1000GEA$MAF == 0)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.2389546
## two counts
(sum(mafin1000GEA$MAF <= 2/498/2)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3260734

tmp <- marginEA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginAA$BETA) ) # 0.9601227
cor(tmp,marginAA$BETA) # 0.9561637

df.AA <- data.frame(Beta_EA=tmp,Beta_AA=marginAA$BETA,ID=marginAA$ID,stringsAsFactors = F)

df.AA[(df.AA$Beta_AA < -0.5) & (df.AA$Beta_EA > 0.15),] # 893
df.AA[(df.AA$Beta_AA > 0.5) & (df.AA$Beta_EA < -0.2),] # 59
df.AA[(df.AA$Beta_AA > 1) & (df.AA$Beta_EA < 1),] # 710
df.AA[(df.AA$Beta_AA > 0) & (df.AA$Beta_AA < 0.5) & (df.AA$Beta_EA > 1),] # 168

df.AA_highlight <- df.AA[c(893,59,710,168), ]

p.AA <- ggplot(data = df.AA, aes(x = Beta_AA, y = Beta_EA)) + 
  geom_point(size=0.5,col="#238b45") +
  geom_abline(intercept = 0, slope = 1, col="red") +
  theme(axis.line = element_line(color="black", size = 0.2)) +
  ylim(-2,2)+xlim(-2,2)+
  mdthemes::md_theme_classic() +
  labs(x = "Effect size (AA)", 
       y = "Effect size (EA)",
       title="Common sentinel SNP of *cis*-pQTLs in AA"
       ) +
  My_Theme +
  geom_point(data=df.AA_highlight, aes(x = Beta_AA, y = Beta_EA), size=0.5, col="darkorange")


###############################################################
###############################################################
###############################################################


p <- cowplot::plot_grid(p.EA, p.AA, ncol=2)

ggsave(filename="sp1.pdf",
       plot=p, device="pdf",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/",
       width=160, height=80, units="mm", dpi=320)




