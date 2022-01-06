
## Suppl Fig 5

library(readr)
library(ggplot2)
library(ggpubr)

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

# V7

## AA Liver
df.hsqAALiver <- read_excel("ExtendedDataFig5.xlsx", sheet = "5a-1")
df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p11 <- ggplot(data = df.hsqAALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqAALiver$gene))," overlapping genes\nfor T in liver (V7)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## AA Whole_Blood
df.hsqAAWhole_Blood <- read_excel("ExtendedDataFig5.xlsx", sheet = "5a-2")
df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p12 <- ggplot(data = df.hsqAAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqAAWhole_Blood$gene))," overlapping genes\nfor T in whole blood (V7)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))



## EA Liver
df.hsqEALiver <- read_excel("ExtendedDataFig5.xlsx", sheet = "5a-3")
df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p13 <- ggplot(data = df.hsqEALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqEALiver$gene))," overlapping genes\nfor T in liver (V7)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## EA Whole_Blood
df.hsqEAWhole_Blood <- read_excel("ExtendedDataFig5.xlsx", sheet = "5a-4")
df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p14 <- ggplot(data = df.hsqEAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqEAWhole_Blood$gene))," overlapping genes\nfor T in whole blood (V7)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))



###############################################################
###############################################################
###############################################################

# V8

## AA Liver
df.hsqAALiver <- read_excel("ExtendedDataFig5.xlsx", sheet = "5b-1")
df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p21 <- ggplot(data = df.hsqAALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqAALiver$gene))," overlapping genes\nfor T in liver (V8)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## AA Whole_Blood
df.hsqAAWhole_Blood <- read_excel("ExtendedDataFig5.xlsx", sheet = "5b-2")
df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p22 <- ggplot(data = df.hsqAAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqAAWhole_Blood$gene))," overlapping genes\nfor T in whole blood (V8)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))



## EA Liver
df.hsqEALiver <- read_excel("ExtendedDataFig5.xlsx", sheet = "5b-3")
df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p23 <- ggplot(data = df.hsqEALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqEALiver$gene))," overlapping genes\nfor T in liver (V8)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## EA Whole_Blood
df.hsqEAWhole_Blood <- read_excel("ExtendedDataFig5.xlsx", sheet = "5b-4")
df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p24 <- ggplot(data = df.hsqEAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(unique(df.hsqEAWhole_Blood$gene))," overlapping genes\nfor T in whole blood (V8)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))

###############################################################
###############################################################
###############################################################

p <- ggarrange(ggarrange(p11, p12, p13, p14,
                         nrow=1, labels = c("a", NA, NA, NA)
                         ),
               ggarrange(p21, p22, p23, p24,
                         nrow=1, labels = c("b", NA, NA, NA)
                         ),
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.5, 0.5)
)

ggsave(filename=paste0("ExtendedDataFigure5.pdf"), 
       plot=p, device="pdf",
       path="/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Figures/", 
       width=180, height=135, units="mm", dpi=320)


