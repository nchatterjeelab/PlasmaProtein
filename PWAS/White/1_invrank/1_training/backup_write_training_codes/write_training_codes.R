
rm(list=ls())

dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_all_protein')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_all_protein')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_all_protein')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp')

library(readr)
library(stringr)

seqid <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/seqid_autosomal_withSNP.txt")
seqid_exist <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp")
seqid_exist <- seqid_exist[str_detect(seqid_exist, ".wgt.RDat")]
seqid_exist <- unlist(lapply(str_split(seqid_exist,".wgt.RDat"), FUN=function (x){x[1]}))
seqid <- setdiff(seqid, seqid_exist)

a <- "#!/usr/bin/env bash
#$ -N submitjobs
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e

module load old_conda_R/3.6

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights_plinkthreads.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq_remove_ambiguous_snp/", seqid[i], " \\
--verbose 2 \\
--hsq_set 0.01 \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp_all_protein/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--models enet \\
--save_hsq TRUE

")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles_all_protein/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles_all_protein/", seqid[i],".sh
")

  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs_all_protein/submitjobs.sh')


####################################
####################################


rm(list=ls())

dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp')
dir.create('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp')

library(readr)
library(stringr)

seqid <- readLines("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/seqid_autosomal_withSNP.txt")

tmp1 <- c("SeqId_10513_13","SeqId_10817_26","SeqId_11837_7","SeqId_11872_9","SeqId_13933_276","SeqId_14112_40","SeqId_14314_6",
"SeqId_15321_8","SeqId_15322_35","SeqId_15386_7","SeqId_15449_33","SeqId_15591_28","SeqId_16308_14","SeqId_16803_4",
"SeqId_16845_15","SeqId_16851_50","SeqId_16854_17","SeqId_16865_62","SeqId_16875_13","SeqId_16892_23","SeqId_16908_5",
"SeqId_16926_44","SeqId_16932_5","SeqId_17152_10","SeqId_17154_2","SeqId_17163_117","SeqId_17166_4","SeqId_17175_5",
"SeqId_17195_43","SeqId_17204_17","SeqId_17209_27","SeqId_17210_2","SeqId_17231_1","SeqId_17319_1","SeqId_17320_19",
"SeqId_17326_44","SeqId_17336_54","SeqId_17337_1","SeqId_17347_80","SeqId_17348_5","SeqId_17355_56","SeqId_17357_33",
"SeqId_17370_186","SeqId_17372_5","SeqId_17384_110","SeqId_17400_71","SeqId_17403_14","SeqId_17404_5","SeqId_17408_2",
"SeqId_17435_43","SeqId_17466_72","SeqId_17474_106","SeqId_17475_18","SeqId_17505_125","SeqId_17510_7","SeqId_17511_10",
"SeqId_17512_2","SeqId_17516_7","SeqId_17683_2","SeqId_17710_40","SeqId_17711_13","SeqId_17721_82","SeqId_17724_3",
"SeqId_17725_37","SeqId_17728_61","SeqId_17729_20","SeqId_17735_130","SeqId_17742_2","SeqId_17746_77","SeqId_17751_68",
"SeqId_17758_79","SeqId_17760_128","SeqId_17761_2","SeqId_17770_42","SeqId_17771_35","SeqId_17776_15","SeqId_17777_31",
"SeqId_17781_191","SeqId_17785_11","SeqId_17787_1","SeqId_17789_1","SeqId_17791_25","SeqId_17792_158","SeqId_17793_4",
"SeqId_17794_6","SeqId_17795_176","SeqId_17797_1","SeqId_17800_33","SeqId_17805_35","SeqId_17811_78","SeqId_17813_21",
"SeqId_17818_22","SeqId_17822_57","SeqId_17823_40","SeqId_17826_341","SeqId_17827_53","SeqId_17828_3","SeqId_17835_28",
"SeqId_17849_6","SeqId_17852_5","SeqId_17855_28","SeqId_17857_6","SeqId_18156_7","SeqId_18158_45","SeqId_18160_2",
"SeqId_18165_181","SeqId_18170_46","SeqId_18172_71","SeqId_18177_49","SeqId_18178_13","SeqId_18181_2","SeqId_18182_24",
"SeqId_18184_28","SeqId_18187_16","SeqId_18188_12","SeqId_18189_12","SeqId_18190_15","SeqId_18191_23","SeqId_18197_97",
"SeqId_18198_51","SeqId_18202_22","SeqId_18204_1","SeqId_18206_18","SeqId_18208_3","SeqId_18210_12","SeqId_18213_30",
"SeqId_18226_148","SeqId_18227_3","SeqId_18228_30","SeqId_18231_147","SeqId_18232_42","SeqId_18237_29","SeqId_18240_6",
"SeqId_18241_18","SeqId_18242_8","SeqId_18255_6","SeqId_18257_64","SeqId_18259_15","SeqId_18261_34","SeqId_18264_12",
"SeqId_18265_18","SeqId_18267_74","SeqId_18271_43","SeqId_18275_5","SeqId_18276_34","SeqId_18277_28","SeqId_18282_1",
"SeqId_18284_77","SeqId_18285_6","SeqId_18286_3","SeqId_18290_6","SeqId_18291_8","SeqId_18294_26","SeqId_18295_102",
"SeqId_18297_8","SeqId_18301_10","SeqId_18302_204","SeqId_18303_39","SeqId_18304_19","SeqId_18306_1","SeqId_18307_71",
"SeqId_18308_30","SeqId_18309_18","SeqId_18310_26","SeqId_18312_68","SeqId_18313_4","SeqId_18315_38","SeqId_18316_75",
"SeqId_18317_111","SeqId_18318_98","SeqId_18319_7","SeqId_18321_38","SeqId_18322_15","SeqId_18323_39","SeqId_18324_61",
"SeqId_18327_6","SeqId_18328_36","SeqId_18330_7","SeqId_18332_17","SeqId_18338_26","SeqId_18339_207","SeqId_18340_2",
"SeqId_18342_2","SeqId_18343_10","SeqId_18348_89","SeqId_18373_13","SeqId_18375_28","SeqId_18376_19","SeqId_18380_78",
"SeqId_18381_16","SeqId_18382_109","SeqId_18386_36","SeqId_18387_7","SeqId_18389_11","SeqId_18392_19","SeqId_18395_5","SeqId_18396_10")

tmp2 <- seqid[(which(seqid=="SeqId_18396_10")+1):length(seqid)]

seqid <- c(tmp1, tmp2)

a <- "#!/usr/bin/env bash
#$ -N submitjobs
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e

module load old_conda_R/3.6

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights_plinkthreads.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq_remove_ambiguous_snp/", seqid[i], " \\
--verbose 2 \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--save_hsq TRUE \\
--models enet


")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles/", seqid[i],".sh
")

  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs/submitjobs.sh')


####################################
####################################

seqid <- c("SeqId_18375_28","SeqId_18380_78","SeqId_17474_106","SeqId_17827_53","SeqId_18261_34",
"SeqId_18324_61","SeqId_18434_141","SeqId_18875_125","SeqId_18918_86","SeqId_19135_5",
"SeqId_19197_95","SeqId_19265_9","SeqId_19347_37","SeqId_19545_145","SeqId_2190_55",
"SeqId_2212_69","SeqId_2247_20","SeqId_2278_61","SeqId_2381_52","SeqId_2421_7",
"SeqId_2431_17","SeqId_2441_2","SeqId_2443_10","SeqId_2475_1","SeqId_2500_2",
"SeqId_2505_49","SeqId_2514_65","SeqId_2515_14","SeqId_2516_57","SeqId_2524_56",
"SeqId_2558_51","SeqId_2567_5","SeqId_2571_12","SeqId_2578_67","SeqId_3516_60",
"SeqId_2642_4","SeqId_2649_77","SeqId_3220_40","SeqId_2666_53","SeqId_2677_1",
"SeqId_2681_23","SeqId_3421_54")

a <- "#!/usr/bin/env bash
#$ -N submitjobs
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

####################
"

for(i in 1:length(seqid)){

b <- paste0("#!/usr/bin/env bash
#$ -N ", seqid[i], "
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
#$ -m e

module load old_conda_R/3.6

Rscript /dcl01/chatterj/data/jzhang2/TWAS/fusion_twas-master/FUSION.compute_weights_plinkthreads.R \\
--PATH_plink /dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \\
--PATH_gcta /dcl01/chatterj/data/jzhang2/TOOLS/gcta_1.93.0beta/gcta64 \\
--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window1M/byseq_remove_ambiguous_snp/", seqid[i], " \\
--verbose 2 \\
--tmp /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/tmp/", seqid[i], " \\
--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_remove_ambiguous_snp/", seqid[i], " \\
--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/peernum/invrankpheno/120/", seqid[i], ".pheno \\
--save_hsq TRUE \\
--models enet


")

  writeLines(b,  paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/shfiles/', seqid[i], '.sh'))


a <- paste0(a,
"
qsub ../shfiles/", seqid[i],".sh
")

  print(i)
}

writeLines(a,  '/dcl01/chatterj/data/jzhang2/pwas/pipeline/codes/GRCh38/PWAS/White/1_invrank/1_training/submitjobs/submitjobs_eqw.sh')

####################################
####################################


a <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs") # coefs_1
b <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun_1") # coefs_2
c <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs_rerun") # coefs_3

length(a) #2389
length(b) #1219
length(c) #51

length(intersect(a,b)) #271
length(intersect(b,c)) #1
length(intersect(a,c)) #30

d <- dir("/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/para1/invrank/coefs") # 3412
