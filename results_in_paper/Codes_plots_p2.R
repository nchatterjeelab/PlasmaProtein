ss = readRDS("~/pQTL/SeqId_6919_3.rds")
ld.aa.inx <- ld.aa[which(a1$ss.aa[,2] == 161106),]
ld.aa = ss$ld.aa
ld.aa.inx <- ld.aa[which(a1$ss.aa[,2] == 161106),]
ld.aa.inx <- ld.aa[which(ss$ss.aa[,2] == 161106),]
ld.aa.inx2 <- ld.aa.inx^2
ld.aa.inx2 <- cut(x = ld.aa.inx2,breaks = c(0,0.2,0.4,0.6,0.8,1),include.lowest = T,right = F)
levels(ld.aa.inx2) <- c("lightgrey","deepskyblue2","darkolivegreen4","goldenrod2","firebrick2")
table(ld.aa.inx2)
plot(ss$ss.aa[1:1000,2],(ss$ss.aa[1:1000,5]),col = as.character(ld.aa.inx2[1:1000]),pch = 19,axes = F)
head(ss$ss.aa)
plot(ss$ss.aa[1:1000,2],-log10(ss$ss.aa[1:1000,6]),col = as.character(ld.aa.inx2[1:1000]),pch = 19,axes = F)
ss.aa = ss$ss.aa
ss.ea = ss$ss.ea
dim(ss.ea)
dim(ss.aa)
length(intersect(ss.ea$V2,ss.aa$V2))
ss.aa.1 = ss.aa[match(intersect(ss.ea$V2,ss.aa$V2),ss.aa$V2),]
dim(ss.aa.1)
head(ss.aa.1)
plot(ss.aa.1[,2],-log10(ss.aa.1[,6]),col = as.character(ld.aa.inx2[1:1000]),pch = 19,axes = F)
ss.aa <- cbind(ss.aa,ld.aa.inx2)
head(ss.aa)
ss.aa.1 = ss.aa[match(intersect(ss.ea$V2,ss.aa$V2),ss.aa$V2),]
ld.ea = ss$ld.ea
ld.ea.inx <- ld.ea[which(ss$ss.ea[,2] == 161106),]
ld.ea.inx2 <- ld.ea.inx^2
ld.ea.inx2 <- cut(x = ld.ea.inx2,breaks = c(0,0.2,0.4,0.6,0.8,1),include.lowest = T,right = F)
levels(ld.ea.inx2) <- c("lightgrey","deepskyblue2","darkolivegreen4","goldenrod2","firebrick2")
ss.ea <- cbind(ss.ea,ld.ea.inx2)



mh1 <- ggplot(ss.ea[1:1000,],aes(x=pos,y=(p),color = ss.ea$cols[1:1000])) + geom_point(alpha = 1,size = 1.8,aes(color = ss.ea$cols[1:1000])) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=10),axis.title.y = element_text(size = 8),legend.position = c(0.75, 0.65),legend.text = element_text(size = 6),legend.title = element_text(size = 7),
        legend.background = element_rect(fill = "white")) +
  scale_x_continuous(name = "",breaks = c(ss.ea$pos[1],ss.ea$pos[sentinel.ea],ss.ea$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.ea$pos[sentinel.ea],ss.ea$pos[1000])/1000000,3))) +
  scale_y_continuous(name = expression("-log"[10]*"(p-value)"),limits = c(0,450),breaks = c(0,200,400)) + 
  annotate("point", x = ss.ea$pos[314], y = ss.ea$p[314], colour = "darkmagenta",size = 3.5,shape = 18) + 
  scale_color_manual(name = expression("r"^2), labels = c("0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1"), values = c("lightgrey","deepskyblue2","darkolivegreen4","goldenrod2","firebrick2")) + 
  annotate(x=ss.ea$pos[100], y=400, geom = "text",label="EA",col = "black",size = 8)

mh2 <- ggplot(ss.aa.1[1:1000,],aes(x=pos,y=-log10(p))) + geom_point(alpha = 1,color = ss.aa.1$cols[1:1000],size = 1.8) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=10),axis.title.y = element_text(size = 8)) +
  scale_x_continuous(name = paste("Position on chromosome:",ss.ea[1,1],"(Mb)"), breaks = c(ss.aa.1$pos[1],ss.aa.1$pos[sentinel.aa],ss.aa.1$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.aa.1$pos[sentinel.aa],ss.ea$pos[1000])/1000000,3))) +
  scale_y_continuous(name = expression("-log"[10]*"(p-value)"),limits = c(0,130),breaks = c(0,100),) + 
  annotate("point", x = ss.aa.1$pos[262], y = -log10(ss.aa.1$p[262]), colour = "darkmagenta",size = 3.5,shape = 18) + 
  annotate(x=ss.aa.1$pos[100], y=100, geom = "text",label="AA",col = "black",size = 8)
# mh2 <- mh2+My_Theme


sent.ea <- which(ss.ea.1$pos == 161106)
mh3 <- ggplot(ss.ea.1[1:1000,],aes(x=pos,y=(pip))) + geom_point(alpha = 1,color = ss.ea.1$cols[1:1000],size = 1.8) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=10),axis.title.y = element_text(size = 8)) +
  scale_x_continuous(name = "",breaks = c(ss.ea.1$pos[1],ss.ea.1$pos[sent.ea],ss.ea.1$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.ea.1$pos[sent.ea],ss.ea$pos[1000])/1000000,3))) +
  scale_y_continuous(name = "PIP",limits = c(0,1),breaks = c(0,0.5,1)) + 
  annotate("point", x = ss.ea.1$pos[262], y = ss.ea.1$pip[262], colour = "darkmagenta",size = 2.5,shape = 18) + 
  annotate(x=ss.ea.1$pos[100], y=0.9, geom = "text",label="EA",col = "black",size = 8)
# mh3 <- mh3+My_Theme


# sentinel.aa <- which(ss.aa$pos == 161106)
mh4 <- ggplot(ss.aa.1[1:1000,],aes(x=pos,y=(pip))) + geom_point(alpha = 1,color = ss.aa.1$cols[1:1000],size = 1.8) +
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.2),axis.text = element_text(size=10),axis.title.y = element_text(size = 8)) +
  scale_x_continuous(name = paste("Position on chromosome:",ss.ea[1,1],"(Mb)"), breaks = c(ss.aa.1$pos[1],ss.aa.1$pos[sentinel.aa],ss.aa.1$pos[1000]),labels = as.character(round(c(ss.ea$pos[1],ss.aa.1$pos[sentinel.aa],ss.ea$pos[1000])/1000000,3))) +
  scale_y_continuous(name = "PIP",limits = c(0,1),breaks = c(0,0.5,1)) + 
  annotate("point", x = ss.aa.1$pos[262], y = ss.aa.1$pip[262], colour = "darkmagenta",size = 2.5,shape = 18) + 
  annotate(x=ss.aa.1$pos[100], y=0.8, geom = "text",label="AA",col = "black",size = 8)
# mh4 <- mh4+My_Theme


ll1 = readRDS("~/pQTL/cs_plot.rds")
as3 = cbind(ll1[,1]-6.5,ll1[,2]+2)
as3 = data.frame(c(as3[,1],as3[,2]),rep(c("EA","AA"),each = dim(as3)[1]))
bx1 <- ggplot(as3, aes(x=as3[,2], y=as3[,1],fill = as3[,2])) +
  geom_boxplot(alpha = 0.8,notch = T,notchwidth = 0.5,outlier.shape = NA)+     scale_fill_manual(values=c("#238b45","#2171b5")) +
  coord_cartesian(ylim = c(0,50)) + theme(plot.title = element_text(hjust = 0.5,size = 10),
                                          panel.background = element_blank(),
                                          axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size = 12),
                                          legend.position="none")+
  labs(title = "Variants in credible sets", x=NULL,y=NULL)
# bx1 <- bx1+My_Theme

ll2 = readRDS("~/pQTL/susie.rds")
as4 = data.frame(c(ll2[[1]][,1],ll2[[2]][,1]),rep(c("EA","AA"),each = 1447))
bx2 <- ggplot(as4, aes(x=as4[,2], y=as4[,1],fill = as4[,2])) +
  geom_boxplot(alpha = 0.8,notch = T,notchwidth = 0.5,outlier.shape = NA)+     scale_fill_manual(values=c("#238b45","#2171b5")) +
  coord_cartesian(ylim = c(0,10)) + theme(plot.title = element_text(hjust = 0.5,size = 10),
                                          panel.background = element_blank(),
                                          axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size = 12),
                                          legend.position="none")+
  labs(title = "Number of components", x=NULL,y=NULL)
# bx2 <- bx2+My_Theme

p1 <- ggarrange(ggarrange(bx1,bx2,ncol = 2,widths = c(0.5,0.5),labels = c("a","b")),ggarrange(mh1,mh3,ncol = 2,widths = c(0.5,0.5),labels = c("c","d")),ggarrange(mh2,mh4,ncol = 2,widths = c(0.5,0.5),labels = c("e","f")),nrow = 3)
