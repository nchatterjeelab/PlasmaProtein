library(ggplot2)
library(stringr)
library(readr)
library(dplyr)


res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/correlations_v7.txt")
gene <- unique(res$gene)
tissue <- unique(res$tissue)
cor <- matrix(nrow=length(gene), ncol=length(tissue))
for(i in 1:length(gene)){
  for (j  in 1:length(tissue)) {
    m <- (res$tissue==tissue[j]) & (res$gene==gene[i])
    if(sum(m) == 0){cor[i,j] <- NA; next}
    cor[i,j] <- res$cor[m]
  }
}
colnames(cor) <- tissue
df <- data.frame(gene=gene, cor)
write_tsv(df, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/7_correlation_tab_v7.txt")

res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/correlations_v8.ALL.txt")
gene <- unique(res$gene)
tissue <- unique(res$tissue)
cor <- matrix(nrow=length(gene), ncol=length(tissue))
for(i in 1:length(gene)){
  for (j  in 1:length(tissue)) {
    m <- (res$tissue==tissue[j]) & (res$gene==gene[i])
    if(sum(m) == 0){cor[i,j] <- NA; next}
    cor[i,j] <- res$cor[m]
  }
}
colnames(cor) <- tissue
df <- data.frame(gene=gene, cor)
write_tsv(df, "/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Tables/7_correlation_tab_v8.txt")




a <- unlist(lapply(strsplit(res$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
res$tissue <- gtex.colors$V1[match(a, b)]
tmp <- res %>% group_by(tissue) %>% summarise(med = median(cor, na.rm=T))
tmp <- tmp[order(tmp$med,decreasing=T),]
res$tissue = factor(res$tissue, levels=tmp$tissue)

myColors <- gtex.colors$V2
names(myColors) <- gtex.colors$V1
colScale <- scale_fill_manual(name = "gtex.colors", values = myColors)

p <- ggplot(data = res, aes(x = tissue, fill=tissue)) +
  geom_boxplot(alpha=0.8, notch = TRUE, notchwidth = 0.5, aes(y=cor)) + 
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none",
        axis.title.y = element_text(hjust=1))+
  coord_cartesian(ylim = c(-0.25,1))+
  colScale +
  labs(x = "GTEx V7 tissue", 
       y = "Correlation between cis-regulated gene       \nexpression and plasma protein SOMAmers      ",
       title=NULL)
# ignore -- (Color represents for the proportion of genes whose expression levels and plasma protein levels are both significant cis-heritable. The lighter the color, the higher the proportion.)

### Put them together
library(ggpubr)

ggsave(filename="p2.png", 
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/", 
       width=12, height=4.6, units="in", dpi=500)


