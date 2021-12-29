
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

annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/prot.anno_autosomal.txt")

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")

tmp <- summary.AA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.AA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]
tmp <- summary.EA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.EA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]

annota.Liver <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/Liver.P01.pos")
annota.Whole_Blood <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/Whole_Blood.P01.pos")

summary.Liver <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Liver.rds")
summary.Whole_Blood <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Whole_Blood.rds")

summary.Liver$gene <- annota.Liver$ID[match(paste0("Liver/",summary.Liver$RDat), annota.Liver$WGT)]
summary.Whole_Blood$gene <- annota.Whole_Blood$ID[match(paste0("Whole_Blood/",summary.Whole_Blood$RDat), annota.Whole_Blood$WGT)]

geneEALiver <- unique(intersect(summary.EA$gene, summary.Liver$gene))
geneEAWhole_Blood <- unique(intersect(summary.EA$gene, summary.Whole_Blood$gene))

geneAALiver <- unique(intersect(summary.AA$gene, summary.Liver$gene))
geneAAWhole_Blood <- unique(intersect(summary.AA$gene, summary.Whole_Blood$gene))

## AA Liver
df.hsqAALiver <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneAALiver]), 
                            group=c(rep("P in AA", sum(summary.AA$gene %in% geneAALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneAALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneAALiver)) ),
                            stringsAsFactors = F)

df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p11 <- ggplot(data = df.hsqAALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneAALiver)," overlapping genes\nfor T in liver (V7)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## AA Whole_Blood
df.hsqAAWhole_Blood <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneAAWhole_Blood]), 
                                  group=c(rep("P in AA", sum(summary.AA$gene %in% geneAAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p12 <- ggplot(data = df.hsqAAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneAAWhole_Blood)," overlapping genes\nfor T in whole blood (V7)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))



## EA Liver
df.hsqEALiver <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneEALiver]), 
                            group=c(rep("P in EA", sum(summary.EA$gene %in% geneEALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneEALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneEALiver)) ),
                            stringsAsFactors = F)

df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p13 <- ggplot(data = df.hsqEALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneEALiver)," overlapping genes\nfor T in liver (V7)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## EA Whole_Blood
df.hsqEAWhole_Blood <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneEAWhole_Blood]), 
                                  group=c(rep("P in EA", sum(summary.EA$gene %in% geneEAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p14 <- ggplot(data = df.hsqEAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneEAWhole_Blood)," overlapping genes\nfor T in whole blood (V7)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))



###############################################################
###############################################################
###############################################################

# V8

annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/prot.anno_autosomal.txt")

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")

tmp <- summary.AA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.AA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]
tmp <- summary.EA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.EA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]

annota.Liver <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/GTEXv8.Liver.geneid.pos")
annota.Whole_Blood <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ge_pos/GTEXv8.Whole_Blood.geneid.pos")

summary.Liver <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v8.ALL/summary.Liver.rds")
summary.Whole_Blood <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v8.ALL/summary.Whole_Blood.rds")

summary.Liver$gene <- annota.Liver$ID[match(summary.Liver$RDat, annota.Liver$WGT)]
summary.Whole_Blood$gene <- annota.Whole_Blood$ID[match(summary.Whole_Blood$RDat, annota.Whole_Blood$WGT)]

geneEALiver <- unique(intersect(summary.EA$gene, summary.Liver$gene))
geneEAWhole_Blood <- unique(intersect(summary.EA$gene, summary.Whole_Blood$gene))

geneAALiver <- unique(intersect(summary.AA$gene, summary.Liver$gene))
geneAAWhole_Blood <- unique(intersect(summary.AA$gene, summary.Whole_Blood$gene))

## AA Liver
df.hsqAALiver <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneAALiver]), 
                            group=c(rep("P in AA", sum(summary.AA$gene %in% geneAALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneAALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneAALiver)) ),
                            stringsAsFactors = F)

df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p21 <- ggplot(data = df.hsqAALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneAALiver)," overlapping genes\nfor T in liver (V8)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## AA Whole_Blood
df.hsqAAWhole_Blood <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneAAWhole_Blood]), 
                                  group=c(rep("P in AA", sum(summary.AA$gene %in% geneAAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p22 <- ggplot(data = df.hsqAAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneAAWhole_Blood)," overlapping genes\nfor T in whole blood (V8)\nand P in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))



## EA Liver
df.hsqEALiver <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneEALiver]), 
                            group=c(rep("P in EA", sum(summary.EA$gene %in% geneEALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneEALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneEALiver)) ),
                            stringsAsFactors = F)

df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p23 <- ggplot(data = df.hsqEALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneEALiver)," overlapping genes\nfor T in liver (V8)\nand P in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  My_Theme +
  scale_fill_manual(values=c("#4a1486","#cb181d"))


## EA Whole_Blood
df.hsqEAWhole_Blood <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneEAWhole_Blood]), 
                                  group=c(rep("P in EA", sum(summary.EA$gene %in% geneEAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p24 <- ggplot(data = df.hsqEAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneEAWhole_Blood)," overlapping genes\nfor T in whole blood (V8)\nand P in EA")) +
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
                         nrow=1, labels = c("a", NA, NA, NA)
                         ),
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.5, 0.5)
)

ggsave(filename=paste0("sp5.pdf"), 
       plot=p, device="pdf",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/", 
       width=180, height=135, units="mm", dpi=320)


