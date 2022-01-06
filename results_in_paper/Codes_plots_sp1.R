
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

df <- read_tsv("ExtendedDataFig1.txt")

###############################################################
###############################################################
###############################################################

# EA (panel a)

df.EA <- df[df$eth=="EA",]
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

df.AA <- df[df$eth=="AA",]
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

ggsave(filename="ExtendedDataFigure1.pdf",
       plot=p, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/",
       width=160, height=80, units="mm", dpi=320)




