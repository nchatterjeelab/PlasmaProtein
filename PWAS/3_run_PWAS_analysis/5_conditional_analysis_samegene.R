

########################################################################
## V7

"
mkdir /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v7.samegene/

Rscript /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/scripts/PWAS.conditional.samegene_CI.R \
--PWAS /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/PWAS_CI_hg19.out \
--TWAS /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/TWAS_CI_v7.out/ \
--tissue_list /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V7_tissue_list.txt \
--tissue_n_gene /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V7_n_gene.rds \
--imputed_P /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_Plasma_Protein.txt \
--imputed_T /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_FUSION_v7/ \
--out /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v7.samegene/
"


### merge all tissues (conditional analysis)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V7_tissue_list.txt")

library(dplyr)
library(readr)
library(stringr)
results <- tibble()
for (tissue in tissue_list){

  results <- rbind(results,
                   cbind(tissue=tissue,read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v7.samegene/Table/",tissue,".txt"), col_types=cols())))
  print(tissue)
}
results$PWAS_hit <- factor(results$PWAS_hit, levels=unique(results$PWAS_hit))
results <- results[order(results$PWAS_hit),]

write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v7.samegene/Table/all-cleaned.out"))







########################################################################
## V8

"
mkdir /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v8.samegene/

Rscript /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/scripts/PWAS.conditional.samegene_CI.R \
--PWAS /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/PWAS_CI_hg38.out \
--TWAS /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/TWAS_CI_v8.out/ \
--tissue_list /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_tissue_list.txt \
--tissue_n_gene /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_n_gene.rds \
--imputed_P /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_Plasma_Protein.txt \
--imputed_T /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/1000G_imputed_EA/1000G_imputed_FUSION_v8.ALL/ \
--out /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v8.samegene/
"


### merge all tissues (conditional analysis)

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GTEx_V8_tissue_list.txt")

library(dplyr)
library(readr)
library(stringr)
results <- tibble()
for (tissue in tissue_list){

  results <- rbind(results,
                   cbind(tissue=tissue,read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v8.samegene/Table/",tissue,".txt"), col_types=cols())))
  print(tissue)
}
results$PWAS_hit <- factor(results$PWAS_hit, levels=unique(results$PWAS_hit))
results <- results[order(results$PWAS_hit),]

write_tsv(results, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Urate/ConditionalAnalysis_CI_v8.samegene/Table/all-cleaned.out"))




