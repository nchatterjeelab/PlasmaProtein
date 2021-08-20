
########################################################################
########################################################################
## plot for direct measured proteins
########################################################################
########################################################################


gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)

res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/correlations_direct_measured_prot_v7.txt")

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

p <- ggplot(data = res, aes(x = tissue, fill=tissue)) +
  geom_boxplot(alpha=0.8, notch = TRUE, notchwidth = 0.5, aes(y=cor)) + 
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none",
        axis.title.y = element_text(hjust=1))+
  coord_cartesian(ylim = c(-0.15,0.45))+
  colScale +
  labs(x = "GTEx V7 tissue", 
       y = "Correlation between cis-regulated gene       \nexpression and plasma protein SOMAmers      ",
       title=NULL)


### Put them together

ggsave(filename="sp6.png", 
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/", 
       width=12, height=4.5, units="in", dpi=500)

