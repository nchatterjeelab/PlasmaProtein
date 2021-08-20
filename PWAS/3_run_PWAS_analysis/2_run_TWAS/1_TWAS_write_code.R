
################################################
################################################
## V7

rm(list=ls())

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V7.txt")

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout/TWAS_CI_v7")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7")

b <- paste0("#!/usr/bin/env bash
#$ -N Gout
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (tissue in tissue_list){
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7/",tissue))


a <- paste0("#!/usr/bin/env bash
#$ -N Gout_",tissue,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/scripts/PWAS.assoc_test_CI.R \\
--sumstats /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GWAS/Gout_EA_cleaned_rsid_BETA_SE_TWAS.txt \\
--weights /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/",tissue,".P01.pos \\
--weights_dir /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v7/ \\
--ref_ld_chr /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/LDREF/1000G.EUR. \\
--force_model enet \\
--chr $SGE_TASK_ID \\
--out /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v7/",tissue,"/chr$SGE_TASK_ID.out

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout/TWAS_CI_v7/",tissue,".sh"))

  b <- paste0(b, "
qsub ",tissue,".sh")

}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout/TWAS_CI_v7/ALL.sh"))




################################################
################################################
## V8

rm(list=ls())

tissue_list <- readLines("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/Tissue_list_GTex_V8.txt")

dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout/TWAS_CI_v8")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8")

b <- paste0("#!/usr/bin/env bash
#$ -N Gout
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (tissue in tissue_list){
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8/",tissue))


a <- paste0("#!/usr/bin/env bash
#$ -N Gout_",tissue,"
#$ -cwd
#$ -t 1-22
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/scripts/PWAS.assoc_test_CI.R \\
--sumstats /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/GWAS/Gout_EA_cleaned_rsid_BETA_SE_TWAS.txt \\
--weights /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.",tissue,".geneid.pos \\
--weights_dir /dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/ \\
--ref_ld_chr /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/LDref/EUR/chr \\
--force_model enet \\
--chr $SGE_TASK_ID \\
--out /dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Results/Gout/TWAS_CI_v8/",tissue,"/chr$SGE_TASK_ID.out

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout/TWAS_CI_v8/",tissue,".sh"))

  b <- paste0(b, "
qsub ",tissue,".sh")

}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/submit/Gout/TWAS_CI_v8/ALL.sh"))



