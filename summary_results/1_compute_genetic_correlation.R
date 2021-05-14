library(dplyr)
library(readr)

tissue_list <- readLines("/Users/jnz/Document/JHU/Research/PWAS/Analysis/Pipeline/PWAS/GTex_V7_tissue_list.txt")

annot.p <- suppressMessages(read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/1000G_imputed/Plasma_Protein_EA_hg19.pos")))
imputed.p <- suppressMessages(read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/1000G_imputed/1000G_imputed_Plasma_Protein.txt")))
m <- match(unlist(lapply(annot.p$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)})),
           colnames(imputed.p))
imputed.p <- imputed.p[,m]

gene.p <- annot.p$ID

res <- tibble()
for (tissue in tissue_list) {
  annot.t <- suppressMessages(read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/1000G_imputed/",tissue,".P01.pos")))
  imputed.t <- suppressMessages(read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/Pipeline/PWAS/PWAS_EA/1000G_imputed_EA/1000G_imputed_FUSION/",tissue,".txt")))
  gene.t <- annot.t$ID
  
  gene <- intersect(unique(gene.p), unique(gene.t))
  cor <- numeric()
  for (i in 1:length(gene)){
    imputed.p_tmp <- imputed.p[,annot.p$ID == gene[i]]
    imputed.t_tmp <- imputed.t[,annot.t$ID == gene[i]]
    
    if(ncol(imputed.p_tmp)>1){
      tmp <- imputed.p_tmp[1]
      for (j in 2:ncol(imputed.p_tmp)) {
        tmp <- tmp+imputed.p_tmp[j]
      }
      imputed.p_tmp <- tmp
    }
    
    if(ncol(imputed.t_tmp)>1){
      tmp <- imputed.t_tmp[1]
      for (j in 2:ncol(imputed.t_tmp)) {
        tmp <- tmp+imputed.t_tmp[j]
      }
      imputed.t_tmp <- tmp
    }
    
    cor[i] <- as.numeric(cor(imputed.p_tmp, imputed.t_tmp))
    
  }
  
  res <- rbind(res, 
               data.frame(tissue=tissue, cor= cor, gene=gene))
  
  print(tissue)
  
}

write_tsv(res,"/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/correlations.txt")

res <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/2_CompareTWASmodels/correlations.txt")
summ <- res %>% group_by(tissue) %>% 
    summarise(ave = mean(cor, na.rm=T), 
              med = median(cor, na.rm = T),
              N=length(gene))
summ <- summ[order(summ$med,decreasing=T),]

gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", 
                          sep = '\t', comment.char = '', stringsAsFactors = F)
a <- unlist(lapply(strsplit(summ$tissue, "_|-"), FUN = function(x){paste(x, collapse = "")}))
b <- unlist(lapply(strsplit(gtex.colors$V1, "-|\\(|\\)| "), FUN = function(x){paste(x, collapse = "")}))
summ$tissue <- gtex.colors$V1[match(a, b)]


Nsample <- read.table(text = "Adipose - Subcutaneous	385	8193
Adipose - Visceral (Omentum)	313	6174
Adrenal Gland	175	4521
Artery - Aorta	267	6463
Artery - Coronary	152	3248
Artery - Tibial	388	8224
Brain - Amygdala	88	1837
Brain - Anterior cingulate cortex (BA24)	109	2710
Brain - Caudate (basal ganglia)	144	3661
Brain - Cerebellar Hemisphere	125	4408
Brain - Cerebellum	154	5855
Brain - Cortex	136	4012
Brain - Frontal Cortex (BA9)	118	3144
Brain - Hippocampus	111	2295
Brain - Hypothalamus	108	2315
Brain - Nucleus accumbens (basal ganglia)	130	3240
Brain - Putamen (basal ganglia)	111	2818
Brain - Spinal cord (cervical c-1)	83	2006
Brain - Substantia nigra	80	1604
Breast - Mammary Tissue	251	5044
Cells - EBV-transformed lymphocytes	117	2758
Cells - Transformed fibroblasts	300	7353
Colon - Sigmoid	203	4874
Colon - Transverse	246	5316
Esophagus - Gastroesophageal Junction	213	4888
Esophagus - Mucosa	358	8061
Esophagus - Muscularis	335	7773
Heart - Atrial Appendage	264	5671
Heart - Left Ventricle	272	5082
Liver	153	2914
Lung	383	7776
Minor Salivary Gland	85	1822
Muscle - Skeletal	491	7409
Nerve - Tibial	361	9657
Ovary	122	2809
Pancreas	220	5094
Pituitary	157	4402
Prostate	132	2797
Skin - Not Sun Exposed (Suprapubic)	335	7458
Skin - Sun Exposed (Lower leg)	414	8879
Small Intestine - Terminal Ileum	122	2879
Spleen	146	4497
Stomach	237	4456
Testis	225	9253
Thyroid	399	9826
Uterus	101	2135
Vagina	106	2013
Whole Blood	369	6007", sep="\t")

summ$Nsample <- Nsample$V2[match(summ$tissue,Nsample$V1)]


