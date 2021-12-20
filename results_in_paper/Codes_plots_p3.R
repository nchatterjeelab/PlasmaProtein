library(ggplot2)
library(stringr)
library(readr)
library(dplyr)


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

## Fig 2

summary.AA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.black.rds")
summary.EA <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/summary.white.rds")

soma <- intersect(summary.AA$RDat, summary.EA$RDat)

summary(summary.AA$rsq.enet[summary.AA$RDat %in% soma] / summary.AA$hsq.all[summary.AA$RDat %in% soma])
summary(summary.EA$rsq.enet[summary.EA$RDat %in% soma] / summary.EA$hsq.all[summary.EA$RDat %in% soma])


summary(summary.AA$hsq.all)
summary(summary.EA$hsq.all)
median(summary.AA$rsq.enet/summary.AA$hsq.all)
median(summary.EA$rsq.enet/summary.EA$hsq.all)
mean(summary.AA$rsq.enet/summary.AA$hsq.all)
mean(summary.EA$rsq.enet/summary.EA$hsq.all)

col2 <- colorRampPalette(brewer.pal(8, "Set1"))(4)

#### A

summary.Liver <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Liver.rds")
summary.Whole_Blood <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/h2_summary/RESULTS_v7/summary.Whole_Blood.rds")

summary.Whole_Blood$M <- length(summary.Whole_Blood$RDat)
summary.Liver$M <- length(summary.Liver$RDat)
df.hsq <- data.frame(hsq = c(summary.AA$hsq.all, summary.EA$hsq.all,
                             summary.Liver$hsq.all, summary.Whole_Blood$hsq.all), 
                     group=c(rep("P in AA", summary.AA$M), rep("P in EA", summary.EA$M),
                             rep("T in Liver", summary.Liver$M), rep("T in Whole Blood", summary.Whole_Blood$M)),
                     kind=c(rep("Plasma protein\nSOMAmers (P)", summary.AA$M), rep("Plasma protein\nSOMAmers (P)", summary.EA$M),
                            rep("Gene expr-\nessions (T)", summary.Liver$M), rep("Gene expr-\nessions (T)", summary.Whole_Blood$M)),
                     stringsAsFactors = F)
df.hsq$group <- factor(df.hsq$group, levels = c("T in Liver", "T in Whole Blood", "P in AA", "P in EA"))
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



##### B

library(RColorBrewer)

df.acc <- data.frame(acc = c(summary.AA$rsq.enet/summary.AA$hsq.all, summary.EA$rsq.enet/summary.EA$hsq.all,
                             summary.Liver$rsq.enet/summary.Liver$hsq.all, summary.Whole_Blood$rsq.enet/summary.Whole_Blood$hsq.all,
                             summary.AA$rsq.top1/summary.AA$hsq.all, summary.EA$rsq.top1/summary.EA$hsq.all,
                             summary.Liver$rsq.top1/summary.Liver$hsq.all, summary.Whole_Blood$rsq.top1/summary.Whole_Blood$hsq.all), 
                     group=c(rep("P in AA", summary.AA$M), rep("P in EA", summary.EA$M),
                             rep("T in Liver", summary.Liver$M), rep("T in Whole Blood", summary.Whole_Blood$M),
                             rep("P in AA", summary.AA$M), rep("P in EA", summary.EA$M),
                             rep("T in Liver", summary.Liver$M), rep("T in Whole Blood", summary.Whole_Blood$M)),
                     Model=c(rep("Elastic Net", summary.AA$M+summary.EA$M+summary.Liver$M+summary.Whole_Blood$M),
                             rep("Best cis-SNP", summary.AA$M+summary.EA$M+summary.Liver$M+summary.Whole_Blood$M)), 
                     stringsAsFactors = F)
df.acc$group <- factor(df.acc$group, levels = c("T in Liver", "T in Whole Blood", "P in AA", "P in EA"))
df.acc$Model <- factor(df.acc$Model, levels=c("Best cis-SNP", "Elastic Net"))

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


acc.EA <- summary.EA$rsq.enet/summary.EA$hsq.all; acc.EA <- acc.EA[!is.na(acc.EA)]
acc.AA <- summary.AA$rsq.enet/summary.AA$hsq.all; acc.AA <- acc.AA[!is.na(acc.AA)]

mood.test(acc.EA, acc.AA)

#### C

pred_acc <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/cross_ethnic/EAtoAA_Predicted_Accuracy.txt")
pred_acc$acc <- pred_acc$rsq/pred_acc$hsq
pred_acc$acc.cross <- pred_acc$rsq.cross/pred_acc$hsq
pred_acc$ratio <- pred_acc$rsq.cross/pred_acc$rsq
df.EAtoAA <- data.frame(acc = c(pred_acc$acc, pred_acc$acc.cross),
                        model=c(rep("AA Model", dim(pred_acc)[1]), 
                                rep("EA Model", dim(pred_acc)[1] )))
mean(pred_acc$acc) # 0.7106388
mean(pred_acc$acc.cross) # 0.3093975

pred_acc <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/cross_ethnic/AAtoEA_Predicted_Accuracy.txt")
pred_acc$acc <- pred_acc$rsq/pred_acc$hsq
pred_acc$acc.cross <- pred_acc$rsq.cross/pred_acc$hsq
pred_acc$ratio <- pred_acc$rsq.cross/pred_acc$rsq
df.AAtoEA <- data.frame(acc = c(pred_acc$acc, pred_acc$acc.cross),
                        model=c(rep("EA Model", dim(pred_acc)[1]), 
                                rep("AA Model", dim(pred_acc)[1] )))
mean(pred_acc$acc) # 0.9319751
mean(pred_acc$acc.cross) # 0.5067742

df.CE <- cbind(rbind(df.EAtoAA, df.AAtoEA), 
               race=c(rep("AA",nrow(df.EAtoAA)),rep("EA",nrow(df.AAtoEA))))

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


#### D

# n_Gene <- readRDS("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/n_gene.rds")
# all <- dir("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/")
# all <- all[str_detect(all,"txt")]
# for (i in 1:length(all)) {
#     tmp <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/", all[i]))
#     if(i==1){
#         res <- data.frame(cor=tmp$cor, tissue=all[i], 
#                           geneinp=dim(tmp)[1]/n_Gene[gsub("[.].*", "", all[i])], 
#                           stringsAsFactors = F)
#     }else{
#         res <- rbind(res,
#                      data.frame(cor=tmp$cor, tissue=all[i], 
#                                 geneinp=dim(tmp)[1]/n_Gene[gsub("[.].*", "", all[i])]), 
#                      stringsAsFactors = F)
#     }
# }
# res$tissue <- unlist(strsplit(res$tissue, ".txt"))

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)

res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/correlations_v7.txt")

# a <- unlist(lapply(strsplit(unlist(strsplit(res$tissue, ".txt")), "_|-"), FUN = function(x){paste(x, collapse = "")}))
a <- unlist(lapply(strsplit(res$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
res$tissue <- gtex.colors$V1[match(a, b)]
tmp <- res %>% group_by(tissue) %>% summarise(med = median(cor, na.rm=T))
tmp <- tmp[order(tmp$med,decreasing=T),]
res$tissue = factor(res$tissue, levels=tmp$tissue)

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
# ignore -- (Color represents for the proportion of genes whose expression levels and plasma protein levels are both significant cis-heritable. The lighter the color, the higher the proportion.)

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
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/",
       width=200, height=115, units="mm", dpi=320)


