
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)

tissue <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland",
            "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala",
            "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
            "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex",
            "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus",
            "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
            "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue",
            "Blood_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid",
            "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa",
            "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle",
            "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial",
            "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic",
            "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen",
            "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")


a <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -m e

####################
")

for (i in 1:length(tissue)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", tissue[i],"
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=10G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

tissue=",tissue[i],"

cd /dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/2_CompareTWASmodels/shfiles

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../2_CompareTWASmodels.R ${tissue}.Rout
}
runr \"--args tissue='${tissue}'\"


")

  a <- paste0(a, "
qsub ", tissue[i],".sh
")
  print(tissue[i])

  writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/2_CompareTWASmodels/shfiles/", tissue[i],".sh"))

}
writeLines(a,  paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/2_CompareTWASmodels/shfiles/all.sh"))
