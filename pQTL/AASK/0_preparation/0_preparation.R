
library(haven)
#library(data.table)
library(readr)

# Cohort
aask_baseline <- read_dta('/dcl02/leased/kidney/AASK/static/Cohort/baseline_new.dta')
aask_baseline <- aask_baseline[c('pid','female','age_r')]
aask_baseline$gender <- ifelse(aask_baseline$female == 1, 2, 1)
aask_baseline <- aask_baseline[,c("pid","age_r","gender")]

# Genetic PCs
PCs <- bigreadr::fread2(file="/dcl02/leased/kidney/AASK/static/Genetics/metadata/genetic_pcs.txt")
PCs$gene_id <- paste0("0_",PCs$gene_id)
PCs <- PCs[,c("gene_id","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]

# CRIC ID Mapping
id_mapping <- read.table("/dcl02/leased/kidney/AASK/static/Genetics/metadata/id_linkage.txt",header = T, stringsAsFactors = F,sep = ",")
id_mapping$gene_id <- paste0("0_",id_mapping$gene_id)

merged <- merge(aask_baseline,id_mapping,by.x="pid",by.y="pid",all.x=TRUE, all.y=TRUE)
merged <- merge(merged, PCs, by.x = "gene_id", by.y = "gene_id", all.x=TRUE)
merged <- merged[!(is.na(merged$PC1)),]
cov <- merged

soma <- bigreadr::fread2("/dcl02/leased/kidney/AASK/static/Proteomics/soma_ver4_log2_protein_SMP_705.txt")

#seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/seqid_autosomal.txt')
#seqid <- seqid[seqid %in% colnames(soma)]
#writeLines(seqid, '/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt')
seqid <- readLines('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/seqid_autosomal_overlapwithARIC.txt')

prot <- soma[,c("pid",seqid)]
merged <- merge(id_mapping, prot,by.x="pid",by.y="pid", all.y=TRUE)
merged <- merged[!is.na(merged$gene_id),]
prot <- merged

gene_id <- intersect(cov$gene_id, prot$gene_id)

prot <- prot[match(gene_id,prot$gene_id),]
cov <- cov[match(gene_id, cov$gene_id),]

write_tsv(prot, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/prot.txt'))
write_tsv(cov, paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/aask.cov'))



