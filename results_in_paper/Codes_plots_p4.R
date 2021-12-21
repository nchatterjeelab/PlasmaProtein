
## Fig 4

###############################################################
###############################################################
###############################################################

# Urate (panel a)

library(readr)
library(ggplot2)
library(ggpubr)
library(latex2exp)


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

disease <- "Urate"

#####################
## data preparation

dat.pwas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_",disease,"_CI_hg19.txt"))
p.pwas <- 0.05/nrow(dat.pwas)
dat.pwas <- dat.pwas[!(is.na(dat.pwas$PWAS.P)),]
dat.pwas <- dat.pwas[!(is.na(dat.pwas$P0)),]
# pos_lookup <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/Plasma_Protein_EA_hg19.pos")
# pos_lookup <- unique(pos_lookup[,c(3,5,6)])
# dat.pwas <- dat.pwas[dat.pwas$ID %in% pos_lookup$ID[!is.na(pos_lookup$P0)],]
# m <- match(dat.pwas$ID,pos_lookup$ID)
# dat.pwas$P0 <- pos_lookup$P0[m]
# dat.pwas$P1 <- pos_lookup$P1[m]
dat.pwas <- dat.pwas
dat.pwas$PWAS.P[dat.pwas$ID=="INHBC"] = 10^(-30)
dat.pwas$ID[dat.pwas$ID=="INHBC"] = "INHBC(7.95e-63)"

dat.pwas <- dat.pwas[,c("ID","CHR","P0","PWAS.P")]
colnames(dat.pwas) <- c("ID","CHR","BP","P")
dat.pwas$tissue <- rep("Plasma",nrow(dat.pwas))

dat.pwas <- dat.pwas[!(is.na(dat.pwas$P)),]
p.pwas <- 0.05/nrow(dat.pwas)


library(dplyr)
tissue_list <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/Pipeline/PWAS/GTEx_V7_tissue_list.txt")
res <- tibble()
for(tissue in tissue_list){
  dat.twas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/TWAS/", disease, "_EA/", tissue, ".out"))
  p.twas <- 0.05/nrow(dat.twas)
  dat.twas <- dat.twas[!(is.na(dat.twas$TWAS.P)),]
  dat.twas$tissue <- rep(tissue,nrow(dat.twas))
  res <- rbind(res,dat.twas)
}

dat.twas <- res[,c("ID","CHR","P0","TWAS.P","tissue")]
colnames(dat.twas) <- c("ID","CHR","BP","P","tissue")

dat.twas <- dat.twas[!(is.na(dat.twas$P)),]
p.twas <- 0.05/nrow(dat.twas)

dat_all <- rbind(dat.pwas, dat.twas)



nCHR <- length(unique(dat_all$CHR))
dat_all$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(dat_all$CHR)){
  nbp[i] <- max(dat_all[dat_all$CHR == i,]$BP)
  dat_all$BPcum[dat_all$CHR == i] <- dat_all$BP[dat_all$CHR == i] + s
  s <- s + nbp[i]
}

axis.set <- dat_all %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

#####################
## PWAS

dat <- dat_all[dat_all$tissue=="Plasma",]

label <- c("INHBB","ITIH1","BTN3A3","INHBA","C11orf68","B3GAT3","INHBC(7.95e-63)","SNUPN","NEO1","FASN")

labels_df.pwas <- data.frame(label=label,
                             logP=-log10(dat$P[match(label,dat$ID)]),
                             BPcum=dat$BPcum[match(label,dat$ID)],
                             CHR=dat$CHR[match(label,dat$ID)])
labels_df.pwas <- labels_df.pwas[order(labels_df.pwas$BPcum),]

manhplot.pwas <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                                 color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.8, size=0.8) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat_all$BPcum),max(dat_all$BPcum))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 32 )) +
  # scale_color_manual(name = "gtex.colors", values = myColors)+
  scale_color_manual(values = rep(c("#4292c6", "#08306b"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.pwas),
             linetype='dashed', col="black", size=0.3) +
  guides(color = F) + 
  labs(x = NULL, 
       title = disease) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    # legend.position = "top",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 5, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 8, vjust = 0.5),
    axis.title = element_text(size=8),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8)
  ) + 
  ggrepel::geom_label_repel(data = labels_df.pwas[1:4,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c(5, 30),
                            nudge_y=0.2,
                            direction = "y",
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.pwas[7,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            nudge_x=0.2*10^8,
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.pwas[c(5,6,8:10),],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c(6, 25),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5)+
  My_Theme


#####################
## TWAS


dat <- dat_all[dat_all$tissue!="Plasma",]

m <- rep(F,nrow(dat))
pos <- labels_df.pwas$BPcum
J <- length(pos)
for (i in 1:nrow(dat)) {
  j = 0
  while ( !(m[i]) & (j<J) ) {
    j=j+1
    if( abs(pos[j] - dat$BPcum[i]) < 500000)
      m[i] <- T
  }
}

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)
a <- unlist(lapply(strsplit(dat$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
dat$tissue <- gtex.colors$V1[match(a, b)]
dat$ID <- paste0(dat$ID,"\n(",dat$tissue,")")

myColors <- gtex.colors$V2
names(myColors) <- gtex.colors$V1

dat$tissue[!m & (dat$CHR %in% ((1:11)*2))] <- "black"
dat$tissue[!m & (dat$CHR %in% ((1:11)*2-1))] <- "grey"
myColors <- c(myColors,"#252525","#969696")
names(myColors)[length(myColors)-1] <- "black"
names(myColors)[length(myColors)] <- "grey"

dat$point_alpha <- 0
dat$point_alpha[m] <- 1

labels_df.twas <- tibble()
for (j in 1:J) {
  tmp <- dat[(dat$BPcum < pos[j] + 500000) & (dat$BPcum > pos[j] - 500000),]
  labels_df.twas <- rbind(labels_df.twas,tmp[which.min(tmp$P),])
}

labels_df.twas <- data.frame(label=labels_df.twas$ID,
                             logP=-log10(labels_df.twas$P),
                             BPcum=labels_df.twas$BPcum,
                             CHR=labels_df.twas$CHR)
labels_df.twas <- labels_df.twas[order(labels_df.twas$BPcum),]

dat <- rbind(dat[!m,],dat[m,])
dat <- dat[dat$P>10^(-195),]

manhplot.twas <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                                 color = as.factor(tissue), size = -log10(P))) +
  geom_point(aes(alpha = point_alpha), size=0.8) + 
  scale_alpha_continuous(range = c(0.3, 1)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat_all$BPcum),max(dat_all$BPcum))) +
  scale_y_reverse(limits=c(195, NA), expand=c(0,0))+
  scale_color_manual(name = "gtex.colors", values = myColors)+
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.twas),
             linetype='dashed', col="black", size=0.3) +
  guides(color = F, alpha = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 8, vjust = 0.5),
    axis.title = element_text(size=8),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8)
  )+
  ggrepel::geom_label_repel(data = labels_df.twas[1,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            nudge_x = -0.5*10^8,
                            nudge_y = -38,
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.twas[2,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c( -300,-60),
                            direction="y",
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.twas[3,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            nudge_y = -50,
                            ylim = c( -205,-50),
                            xlim = c(labels_df.twas$BPcum[2], labels_df.twas$BPcum[3]+10^8),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.twas[4,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            nudge_y = -50,
                            ylim = c( -205,-75),
                            direction="y",
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.twas[5:6,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c( -300,-20),
                            xlim = c(labels_df.twas$BPcum[4], labels_df.twas$BPcum[7]+0.3*10^8),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.8) +
  ggrepel::geom_label_repel(data = labels_df.twas[7:10,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            nudge_y = -30,
                            nudge_x = 1.2*10^8,
                            xlim = c(labels_df.twas$BPcum[7]-0.3*10^8, labels_df.twas$BPcum[10]+5*10^8),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.8) +
  My_Theme


p1 <- cowplot::plot_grid(manhplot.pwas, manhplot.twas, ncol=1, align="v")




###############################################################
###############################################################
###############################################################

# Gout (panel b)

library(readr)
library(ggplot2)

disease <- "Gout"

#####################
## data preparation

dat.pwas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/PWAS/PWAS_",disease,"_CI_hg19.txt"))
p.pwas <- 0.05/nrow(dat.pwas)
dat.pwas <- dat.pwas[!(is.na(dat.pwas$PWAS.P)),]
dat.pwas <- dat.pwas[!(is.na(dat.pwas$P0)),]
# pos_lookup <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/pos/Plasma_Protein_EA_hg19.pos")
# pos_lookup <- unique(pos_lookup[,c(3,5,6)])
# dat.pwas <- dat.pwas[dat.pwas$ID %in% pos_lookup$ID[!is.na(pos_lookup$P0)],]
# m <- match(dat.pwas$ID,pos_lookup$ID)
# dat.pwas$P0 <- pos_lookup$P0[m]
# dat.pwas$P1 <- pos_lookup$P1[m]
dat.pwas <- dat.pwas

dat.pwas <- dat.pwas[,c("ID","CHR","P0","PWAS.P")]
colnames(dat.pwas) <- c("ID","CHR","BP","P")
dat.pwas$tissue <- rep("Plasma",nrow(dat.pwas))

dat.pwas <- dat.pwas[!(is.na(dat.pwas$P)),]
p.pwas <- 0.05/nrow(dat.pwas)


library(dplyr)
tissue_list <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/Pipeline/PWAS/GTEx_V7_tissue_list.txt")
res <- tibble()
for(tissue in tissue_list){
  dat.twas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/TWAS/", disease, "_EA/", tissue, ".out"))
  p.twas <- 0.05/nrow(dat.twas)
  dat.twas <- dat.twas[!(is.na(dat.twas$TWAS.P)),]
  dat.twas$tissue <- rep(tissue,nrow(dat.twas))
  res <- rbind(res,dat.twas)
}

dat.twas <- res[,c("ID","CHR","P0","TWAS.P","tissue")]
colnames(dat.twas) <- c("ID","CHR","BP","P","tissue")

dat.twas <- dat.twas[!(is.na(dat.twas$P)),]
p.twas <- 0.05/nrow(dat.twas)

dat_all <- rbind(dat.pwas, dat.twas)



nCHR <- length(unique(dat_all$CHR))
dat_all$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(dat_all$CHR)){
  nbp[i] <- max(dat_all[dat_all$CHR == i,]$BP)
  dat_all$BPcum[dat_all$CHR == i] <- dat_all$BP[dat_all$CHR == i] + s
  s <- s + nbp[i]
}

axis.set <- dat_all %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)


#####################
## PWAS

dat <- dat_all[dat_all$tissue=="Plasma",]

label <- c("IL1RN","BTN3A3","INHBC")

labels_df.pwas <- data.frame(label=label,
                             logP=-log10(dat$P[match(label,dat$ID)]),
                             BPcum=dat$BPcum[match(label,dat$ID)],
                             CHR=dat$CHR[match(label,dat$ID)])
labels_df.pwas <- labels_df.pwas[order(labels_df.pwas$BPcum),]
manhplot.pwas <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                                 color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.8, size=0.8) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat_all$BPcum),max(dat_all$BPcum))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 24)) +
  scale_color_manual(values = rep(c("#4292c6", "#08306b"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  # geom_segment(aes(x=min(dat_all$BPcum), xend=max(dat_all$BPcum),
  #                  y=-log10(p.pwas), yend=-log10(p.pwas)),
  #              linetype='dashed', col="blue", size=0.3) +
  geom_hline(yintercept = -log10(p.pwas),
             linetype='dashed', col="black", size=0.3) +
  guides(color = F) + 
  labs(x = NULL,
       title = disease) +
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 5, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 8, vjust = 0.5),
    axis.title = element_text(size=8),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8)
  ) + 
  ggrepel::geom_label_repel(data = labels_df.pwas[1:2,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c(min(labels_df.pwas$logP), 30),
                            direction = "y", nudge_y=1,
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  ggrepel::geom_label_repel(data = labels_df.pwas[3,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c(5, 28),
                            direction = "y", 
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) +
  My_Theme


#####################
## TWAS


dat <- dat_all[dat_all$tissue!="Plasma",]
ylim <- abs(floor(log10(min(dat$P)))) + 2 

m <- rep(F,nrow(dat))
pos <- labels_df.pwas$BPcum
J <- length(pos)
for (i in 1:nrow(dat)) {
  j = 0
  while ( !(m[i]) & (j<J) ) {
    j=j+1
    if( abs(pos[j] - dat$BPcum[i]) < 500000)
      m[i] <- T
  }
}

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)
a <- unlist(lapply(strsplit(dat$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
dat$tissue <- gtex.colors$V1[match(a, b)]
dat$ID <- paste0(dat$ID,"\n(",dat$tissue,")")

myColors <- gtex.colors$V2
names(myColors) <- gtex.colors$V1

dat$tissue[!m & (dat$CHR %in% ((1:11)*2))] <- "black"
dat$tissue[!m & (dat$CHR %in% ((1:11)*2-1))] <- "grey"
myColors <- c(myColors,"#252525","#969696")
names(myColors)[length(myColors)-1] <- "black"
names(myColors)[length(myColors)] <- "grey"

dat$point_alpha <- 0
dat$point_alpha[m] <- 1

labels_df.twas <- tibble()
for (j in 1:J) {
  tmp <- dat[(dat$BPcum < pos[j] + 500000) & (dat$BPcum > pos[j] - 500000),]
  labels_df.twas <- rbind(labels_df.twas,tmp[which.min(tmp$P),])
}

labels_df.twas <- data.frame(label=labels_df.twas$ID,
                             logP=-log10(labels_df.twas$P),
                             BPcum=labels_df.twas$BPcum,
                             CHR=labels_df.twas$CHR)
labels_df.twas <- labels_df.twas[order(labels_df.twas$BPcum),]

dat <- rbind(dat[!m,],dat[m,])
manhplot.twas <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                                 color = as.factor(tissue), size = -log10(P))) +
  geom_point(aes(alpha = point_alpha), size=0.8) + 
  scale_alpha_continuous(range = c(0.3, 1)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat_all$BPcum), max(dat_all$BPcum))) +
  scale_y_reverse()+
  scale_color_manual(name = "gtex.colors", values = myColors)+
  scale_size_continuous(range = c(0.5,3)) +
  # geom_segment(aes(x=min(dat_all$BPcum), xend=max(dat_all$BPcum),
  #                  y=-log10(p.twas), yend=-log10(p.twas)),
  #              linetype='dashed', col="blue", size=0.3) +
  geom_hline(yintercept = -log10(p.twas), limits = c(min(dat_all$BPcum), max(dat_all$BPcum)),
             linetype='dashed', col="black", size=0.3) +
  guides(color = F, alpha = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 8, vjust = 0.5),
    axis.title = element_text(size=8),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8)
  )+
  ggrepel::geom_label_repel(data = labels_df.twas,
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col=c("#737373","black","black"),
                            size = 2, segment.size = 0.2,
                            point.padding = 0.3, 
                            direction = "y",
                            ylim = c( -130, -20),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5)+
  My_Theme


p2 <- cowplot::plot_grid(manhplot.pwas, manhplot.twas, ncol=1, align="v")


###############################################################
###############################################################
###############################################################

## TWAS tissue color legends

tmp <- ggplot(dat[!(dat$tissue %in% c("black","grey")), ], aes(x = BPcum, y = -log10(P), 
                                                             color = as.factor(tissue))) +
  geom_point() + 
  scale_color_manual(name = "GTEx tissues in TWAS", values = myColors)+
  theme_minimal() +
  theme(
    legend.key.size = unit(2, "mm"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+ 
  My_Theme+
  guides(color=guide_legend(ncol = 1))

p3 <- as_ggplot(get_legend(tmp))


###############################################################
###############################################################
###############################################################


p <- ggarrange(ggarrange(p1, p2,
                         nrow = 2, labels = c("a", "b"),
                         heights = c(0.6,0.4)),
               p3,
               ncol = 2, 
               labels = c(NA, NA),
               widths = c(0.75,0.25)
               )

ggsave(filename=paste0("p4.pdf"), 
       plot=p, device="pdf",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/", 
       width=200, height=150, units="mm", dpi=320)


# ggsave(filename=paste0("p4.png"), 
#        plot=p, device="png",
#        path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Figures/", 
#        width=200, height=150, units="mm", dpi=320)


