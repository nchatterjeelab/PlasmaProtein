
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)

a <- c("GTEx.Adipose_Subcutaneous.P01.tar.bz2", "GTEx.Adipose_Visceral_Omentum.P01.tar.bz2", "GTEx.Adrenal_Gland.P01.tar.bz2",
      "GTEx.Artery_Aorta.P01.tar.bz2", "GTEx.Artery_Coronary.P01.tar.bz2", "GTEx.Artery_Tibial.P01.tar.bz2",
      "GTEx.Brain_Amygdala.P01.tar.bz2", "GTEx.Brain_Anterior_cingulate_cortex_BA24.P01.tar.bz2",
      "GTEx.Brain_Caudate_basal_ganglia.P01.tar.bz2", "GTEx.Brain_Cerebellar_Hemisphere.P01.tar.bz2",
      "GTEx.Brain_Cerebellum.P01.tar.bz2", "GTEx.Brain_Cortex.P01.tar.bz2", "GTEx.Brain_Frontal_Cortex_BA9.P01.tar.bz2",
      "GTEx.Brain_Hippocampus.P01.tar.bz2", "GTEx.Brain_Hypothalamus.P01.tar.bz2", "GTEx.Brain_Nucleus_accumbens_basal_ganglia.P01.tar.bz2",
      "GTEx.Brain_Putamen_basal_ganglia.P01.tar.bz2", "GTEx.Brain_Spinal_cord_cervical_c-1.P01.tar.bz2", "GTEx.Brain_Substantia_nigra.P01.tar.bz2",
      "GTEx.Breast_Mammary_Tissue.P01.tar.bz2", "GTEx.Cells_EBV-transformed_lymphocytes.P01.tar.bz2", "GTEx.Cells_Transformed_fibroblasts.P01.tar.bz2",
      "GTEx.Colon_Sigmoid.P01.tar.bz2", "GTEx.Colon_Transverse.P01.tar.bz2", "GTEx.Esophagus_Gastroesophageal_Junction.P01.tar.bz2",
      "GTEx.Esophagus_Mucosa.P01.tar.bz2", "GTEx.Esophagus_Muscularis.P01.tar.bz2", "GTEx.Heart_Atrial_Appendage.P01.tar.bz2",
      "GTEx.Heart_Left_Ventricle.P01.tar.bz2", "GTEx.Liver.P01.tar.bz2", "GTEx.Lung.P01.tar.bz2", "GTEx.Minor_Salivary_Gland.P01.tar.bz2",
      "GTEx.Muscle_Skeletal.P01.tar.bz2", "GTEx.Nerve_Tibial.P01.tar.bz2", "GTEx.Ovary.P01.tar.bz2", "GTEx.Pancreas.P01.tar.bz2",
      "GTEx.Pituitary.P01.tar.bz2", "GTEx.Prostate.P01.tar.bz2", "GTEx.Skin_Not_Sun_Exposed_Suprapubic.P01.tar.bz2",
      "GTEx.Skin_Sun_Exposed_Lower_leg.P01.tar.bz2", "GTEx.Small_Intestine_Terminal_Ileum.P01.tar.bz2", "GTEx.Spleen.P01.tar.bz2",
      "GTEx.Stomach.P01.tar.bz2", "GTEx.Testis.P01.tar.bz2", "GTEx.Thyroid.P01.tar.bz2",
      "GTEx.Uterus.P01.tar.bz2", "GTEx.Vagina.P01.tar.bz2", "GTEx.Whole_Blood.P01.tar.bz2")

dir.create("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS/")
dir.create("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS")
for (i in 1:length(a)){

b <- paste0("#!/usr/bin/env bash
#$ -N TWAS
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS

wget http://gusevlab.org/projects/fusion/weights/",a[i],"

tar xjf ",a[i],"

")
 print(i)

 writeLines(b,  paste0("/dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS/", a[i],".sh"))

}

