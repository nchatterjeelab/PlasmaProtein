
#"cp -r /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/ /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum_invrankpheno"
#"cp -r /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum/invrankpheno/ /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/pQTL/peernum_invrankpheno"

################################################################
################################################################
################################################################

rm(list=ls())

library(readr)

ethnic='Black'
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/1_peernum_permutation/")
dir.create("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/1_peernum_permutation/generate_bed/")


for (n_peer in (0:25)*10){

b <- paste0("#!/usr/bin/env bash
#$ -N bed_npeer", n_peer, "_",ethnic,"
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\" ../1_generate_bed.R n_peer_",n_peer,"_",ethnic,".Rout
}
runr \"--args n_peer='",n_peer,"' ethnic='",ethnic,"'\"

")

  writeLines(b,  paste0('/dcs04/nilanjan/data/jzhang2/pwas/pipeline/codes/revision_500Kb/1_peernum_permutation/generate_bed/', n_peer, '_',ethnic,'.sh'))

}