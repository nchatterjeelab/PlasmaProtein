
################################################################
################################################################
################################################################

rm(list=ls())

library(readr)

dir.create("/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/generate_bed/")

#for (n_peer in (0:12)*10){

n_peer <- 120

dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/bed_file"))
dir.create(paste0("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/",n_peer,"/vcf_file"))

b <- paste0("#!/usr/bin/env bash
#$ -N npeer", n_peer, "
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../1_generate_bed.R n_peer_",n_peer,".Rout
}
runr \"--args n_peer='",n_peer,"'\"

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/pQTL/White/White_downsampled_to_black_sample_size/generate_bed/', n_peer, '.sh'))


#}
