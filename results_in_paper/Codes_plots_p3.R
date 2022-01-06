
## Fig 3

library(ggplot2)
library(stringr)
library(readr)
library(dplyr)


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

df.hsq <- read_excel("Fig3.xlsx", sheet = "3a")

p1 <- ggplot(data = df.hsq, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, title=NULL) +
  theme(legend.position="top",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486", "#4a1486", "#cb181d","#cb181d"),
                                   vjust = 0.5, hjust = 0.5, angle = 15))+
  My_Theme+
  scale_fill_manual(values=c("#4a1486","#cb181d"))


###############################################################
###############################################################
###############################################################

##### B

library(RColorBrewer)

df.acc <- read_excel("Fig3.xlsx", sheet = "3b")

p2 <- ggplot(data = df.acc, aes(x = group)) +
  geom_boxplot(alpha=0.8, notch = TRUE, notchwidth = 0.5, aes(y=acc, fill=Model)) + 
  coord_cartesian(ylim = c(0,1.2)) +
  labs(title = NULL, x=NULL,
       y=expression(paste(R^2,"/cis-",h^2))) +
  theme(legend.position="top",
        axis.text.x = element_text(color = c("#4a1486", "#4a1486", "#cb181d","#cb181d"),
                                   vjust = 0.5, hjust = 0.5, angle = 15))+
  My_Theme+
  scale_fill_manual(values=c("#feb24c","#41b6c4"))


###############################################################
###############################################################
###############################################################

#### C

df.CE <- read_excel("Fig3.xlsx", sheet = "3c")

p3 <- ggplot(data = df.CE, aes(x = model)) + 
  geom_boxplot(alpha=0.8, notch = TRUE, notchwidth = 0.5, aes(y=acc, fill=model)) + 
  facet_wrap(~race,  ncol=2)+
  labs(title = NULL, x=NULL,
       y=expression(paste(R^2,"/cis-",h^2))) +
  coord_cartesian(ylim = c(0,1.2))  +
  theme(axis.text.x = element_text(color = c("#238b45", "#2171b5"),
                                   vjust = 0.5, hjust = 0.5, angle = 15),
        legend.position="none") +
  My_Theme+
  scale_fill_manual(values=c("#238b45","#2171b5"))


###############################################################
###############################################################
###############################################################

#### D

res <- read_excel("Fig3.xlsx", sheet = "3d")
gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)

myColors <- gtex.colors$V2
names(myColors) <- gtex.colors$V1
colScale <- scale_fill_manual(name = "gtex.colors", values = myColors)

p4 <- ggplot(data = res, aes(x = tissue, fill=tissue)) +
  geom_boxplot(alpha=0.8, notch = TRUE, notchwidth = 0.5, aes(y=cor)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none",
        axis.title.y = element_text(hjust=1))+
  My_Theme+
  coord_cartesian(ylim = c(-0.25,1))+
  colScale +
  labs(x = "GTEx V7 tissue", 
       y = "Correlation between cis-regulated gene       \nexpression and plasma protein SOMAmers      ",
       title=NULL)


###############################################################
###############################################################
###############################################################

### Put them together
library(ggpubr)
p <- ggarrange(ggarrange(p1, p2,
                         p3,
                         ncol = 3, labels = c("a", "b","c"),
                         widths = c(0.29,0.4,0.31)),
               p4,
               nrow = 2, heights = c(0.4,0.6),
               labels = c(NA,"d"))


ggsave(filename="p3.pdf",
       plot=p, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/",
       width=180, height=105, units="mm", dpi=320)


