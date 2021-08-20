
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

dir.create("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v7")
dir.create("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7")


c <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v7/

")

for (i in 1:length(a)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", a[i],"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7

wget http://gusevlab.org/projects/fusion/weights/",a[i],"

tar -xf ",a[i],"

")
  print(i)


c <- paste0(c, "
qsub ", a[i],".sh
")

  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v7/", a[i],".sh"))

}

writeLines(c,  paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v7/submit.sh"))


################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
dir.create("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR")

a <- c("GTEXv8.EUR.Adipose_Subcutaneous.tar.gz",
"GTEXv8.EUR.Adipose_Visceral_Omentum.tar.gz",
"GTEXv8.EUR.Adrenal_Gland.tar.gz",
"GTEXv8.EUR.Artery_Aorta.tar.gz",
"GTEXv8.EUR.Artery_Coronary.tar.gz",
"GTEXv8.EUR.Artery_Tibial.tar.gz",
"GTEXv8.EUR.Brain_Amygdala.tar.gz",
"GTEXv8.EUR.Brain_Anterior_cingulate_cortex_BA24.tar.gz",
"GTEXv8.EUR.Brain_Caudate_basal_ganglia.tar.gz",
"GTEXv8.EUR.Brain_Cerebellar_Hemisphere.tar.gz",
"GTEXv8.EUR.Brain_Cerebellum.tar.gz",
"GTEXv8.EUR.Brain_Cortex.tar.gz",
"GTEXv8.EUR.Brain_Frontal_Cortex_BA9.tar.gz",
"GTEXv8.EUR.Brain_Hippocampus.tar.gz",
"GTEXv8.EUR.Brain_Hypothalamus.tar.gz",
"GTEXv8.EUR.Brain_Nucleus_accumbens_basal_ganglia.tar.gz",
"GTEXv8.EUR.Brain_Putamen_basal_ganglia.tar.gz",
"GTEXv8.EUR.Brain_Spinal_cord_cervical_c-1.tar.gz",
"GTEXv8.EUR.Brain_Substantia_nigra.tar.gz",
"GTEXv8.EUR.Breast_Mammary_Tissue.tar.gz",
"GTEXv8.EUR.Cells_Cultured_fibroblasts.tar.gz",
"GTEXv8.EUR.Cells_EBV-transformed_lymphocytes.tar.gz",
"GTEXv8.EUR.Colon_Sigmoid.tar.gz",
"GTEXv8.EUR.Colon_Transverse.tar.gz",
"GTEXv8.EUR.Esophagus_Gastroesophageal_Junction.tar.gz",
"GTEXv8.EUR.Esophagus_Mucosa.tar.gz",
"GTEXv8.EUR.Esophagus_Muscularis.tar.gz",
"GTEXv8.EUR.Heart_Atrial_Appendage.tar.gz",
"GTEXv8.EUR.Heart_Left_Ventricle.tar.gz",
"GTEXv8.EUR.Kidney_Cortex.tar.gz",
"GTEXv8.EUR.Liver.tar.gz",
"GTEXv8.EUR.Lung.tar.gz",
"GTEXv8.EUR.Minor_Salivary_Gland.tar.gz",
"GTEXv8.EUR.Muscle_Skeletal.tar.gz",
"GTEXv8.EUR.Nerve_Tibial.tar.gz",
"GTEXv8.EUR.Ovary.tar.gz",
"GTEXv8.EUR.Pancreas.tar.gz",
"GTEXv8.EUR.Pituitary.tar.gz",
"GTEXv8.EUR.Prostate.tar.gz",
"GTEXv8.EUR.Skin_Not_Sun_Exposed_Suprapubic.tar.gz",
"GTEXv8.EUR.Skin_Sun_Exposed_Lower_leg.tar.gz",
"GTEXv8.EUR.Small_Intestine_Terminal_Ileum.tar.gz",
"GTEXv8.EUR.Spleen.tar.gz",
"GTEXv8.EUR.Stomach.tar.gz",
"GTEXv8.EUR.Testis.tar.gz",
"GTEXv8.EUR.Thyroid.tar.gz",
"GTEXv8.EUR.Uterus.tar.gz",
"GTEXv8.EUR.Vagina.tar.gz",
"GTEXv8.EUR.Whole_Blood.tar.gz")

dir.create("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.EUR/")

c <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.EUR/

")

for (i in 1:length(a)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", a[i],"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR

wget https://gtexv8fusion.s3-us-west-1.amazonaws.com/EUR/",a[i],"

tar -xf ",a[i],"

")
  print(i)


c <- paste0(c, "
qsub ", a[i],".sh
")

  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.EUR/", a[i],".sh"))

}

writeLines(c,  paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.EUR/submit.EUR.sh"))


################################################################
################################################################
################################################################

rm(list=ls())

library(readr)
dir.create("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL")

a <- c("GTEXv8.ALL.Adipose_Subcutaneous.tar.gz",
"GTEXv8.ALL.Adipose_Visceral_Omentum.tar.gz",
"GTEXv8.ALL.Adrenal_Gland.tar.gz",
"GTEXv8.ALL.Artery_Aorta.tar.gz",
"GTEXv8.ALL.Artery_Coronary.tar.gz",
"GTEXv8.ALL.Artery_Tibial.tar.gz",
"GTEXv8.ALL.Brain_Amygdala.tar.gz",
"GTEXv8.ALL.Brain_Anterior_cingulate_cortex_BA24.tar.gz",
"GTEXv8.ALL.Brain_Caudate_basal_ganglia.tar.gz",
"GTEXv8.ALL.Brain_Cerebellar_Hemisphere.tar.gz",
"GTEXv8.ALL.Brain_Cerebellum.tar.gz",
"GTEXv8.ALL.Brain_Cortex.tar.gz",
"GTEXv8.ALL.Brain_Frontal_Cortex_BA9.tar.gz",
"GTEXv8.ALL.Brain_Hippocampus.tar.gz",
"GTEXv8.ALL.Brain_Hypothalamus.tar.gz",
"GTEXv8.ALL.Brain_Nucleus_accumbens_basal_ganglia.tar.gz",
"GTEXv8.ALL.Brain_Putamen_basal_ganglia.tar.gz",
"GTEXv8.ALL.Brain_Spinal_cord_cervical_c-1.tar.gz",
"GTEXv8.ALL.Brain_Substantia_nigra.tar.gz",
"GTEXv8.ALL.Breast_Mammary_Tissue.tar.gz",
"GTEXv8.ALL.Cells_Cultured_fibroblasts.tar.gz",
"GTEXv8.ALL.Cells_EBV-transformed_lymphocytes.tar.gz",
"GTEXv8.ALL.Colon_Sigmoid.tar.gz",
"GTEXv8.ALL.Colon_Transverse.tar.gz",
"GTEXv8.ALL.Esophagus_Gastroesophageal_Junction.tar.gz",
"GTEXv8.ALL.Esophagus_Mucosa.tar.gz",
"GTEXv8.ALL.Esophagus_Muscularis.tar.gz",
"GTEXv8.ALL.Heart_Atrial_Appendage.tar.gz",
"GTEXv8.ALL.Heart_Left_Ventricle.tar.gz",
"GTEXv8.ALL.Kidney_Cortex.tar.gz",
"GTEXv8.ALL.Liver.tar.gz",
"GTEXv8.ALL.Lung.tar.gz",
"GTEXv8.ALL.Minor_Salivary_Gland.tar.gz",
"GTEXv8.ALL.Muscle_Skeletal.tar.gz",
"GTEXv8.ALL.Nerve_Tibial.tar.gz",
"GTEXv8.ALL.Ovary.tar.gz",
"GTEXv8.ALL.Pancreas.tar.gz",
"GTEXv8.ALL.Pituitary.tar.gz",
"GTEXv8.ALL.Prostate.tar.gz",
"GTEXv8.ALL.Skin_Not_Sun_Exposed_Suprapubic.tar.gz",
"GTEXv8.ALL.Skin_Sun_Exposed_Lower_leg.tar.gz",
"GTEXv8.ALL.Small_Intestine_Terminal_Ileum.tar.gz",
"GTEXv8.ALL.Spleen.tar.gz",
"GTEXv8.ALL.Stomach.tar.gz",
"GTEXv8.ALL.Testis.tar.gz",
"GTEXv8.ALL.Thyroid.tar.gz",
"GTEXv8.ALL.Uterus.tar.gz",
"GTEXv8.ALL.Vagina.tar.gz",
"GTEXv8.ALL.Whole_Blood.tar.gz")

dir.create("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.ALL/")

c <- paste0("#!/usr/bin/env bash
#$ -N submit
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

cd /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.ALL/

")

for (i in 1:length(a)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", a[i],"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


cd /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL

wget https://gtexv8fusion.s3-us-west-1.amazonaws.com/ALL/",a[i],"

tar -xf ",a[i],"

")
  print(i)


c <- paste0(c, "
qsub ", a[i],".sh
")

  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.ALL/", a[i],".sh"))

}

writeLines(c,  paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/download_WEIGHTS_v8.ALL/submit.ALL.sh"))
