
################################################################
################################################################
################################################################

rm(list=ls())


dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/4_White_matchNblack/all_sample_peers/1_generate_bed/")
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/bed_file"))

for (r in 1:10){

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/bed_file/",r))


  b <- paste0("#!/usr/bin/env bash
#$ -N r", r, "
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
  R CMD BATCH --no-save --no-restore \"$1\" ../1_generate_bed.R r_",r,".Rout
}
runr \"--args r='",r,"'\"

              ")

  writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/4_White_matchNblack/all_sample_peers/1_generate_bed/', r, '.sh'))


}