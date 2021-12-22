
## Suppl Fig 6

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

###############################################################
###############################################################
###############################################################

## plot for direct measured proteins


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
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none",
        axis.title.y = element_text(hjust=1))+
  My_Theme +
  coord_cartesian(ylim = c(-0.15,0.45))+
  colScale +
  labs(x = "GTEx V7 tissue", 
       y = "Correlation between cis-regulated gene       \nexpression and plasma protein SOMAmers      ",
       title=NULL)

###############################################################
###############################################################
###############################################################

ggsave(filename="sp6.pdf", 
       plot=p, device="pdf",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/sp/", 
       width=200, height=75, units="mm", dpi=320)

