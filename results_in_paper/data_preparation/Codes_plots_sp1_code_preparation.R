
## Suppl Fig 1

###############################################################
###############################################################
###############################################################

# EA (panel a)

## EA's pQTLs that are rare in AA
load("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ethnic_specific/results_EA.RData")
(sum(mafin1000GAA$MAF < 0.01)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1334333
(sum(mafin1000GAA$MAF < 0.005)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.1304348
(sum(mafin1000GAA$MAF < 0.002)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09995002
(sum(mafin1000GAA$MAF == 0)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.04197901
## two counts
(sum(mafin1000GAA$MAF <= 2/659/2)) / (nrow(marginAA)+nrow(mafin1000GAA)) # 0.09995002

marginEA$A1 == marginAA$A1
tmp <- marginAA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginEA$BETA) ) # 0.9032445
cor(tmp,marginEA$BETA) # 0.9275709

df.EA <- data.frame(Beta_AA=tmp,Beta_EA=marginEA$BETA,ID=marginEA$ID,stringsAsFactors = F)

df.EA[(df.EA$Beta_AA > 0.2) & (df.EA$Beta_EA < -0.6),] # 865
df.EA[(df.EA$Beta_AA > 0) & (df.EA$Beta_EA < -1),] # 1391
df.EA[(df.EA$Beta_AA < 0) & (df.EA$Beta_EA > 0.5),] # 1342
df.EA[(df.EA$Beta_AA > 0.2) & (df.EA$Beta_AA < 0.3) & (df.EA$Beta_EA > 1),] # 277

df.EA$eth <- "EA"






###############################################################
###############################################################
###############################################################

# AA (panel b)

## AA's pQTLs that are rare in EA
load("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*RData/ethnic_specific/results_AA.RData")
(sum(mafin1000GEA$MAF < 0.01)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3665215
(sum(mafin1000GEA$MAF < 0.005)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3447418
(sum(mafin1000GEA$MAF < 0.002)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3036714
(sum(mafin1000GEA$MAF == 0)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.2389546
## two counts
(sum(mafin1000GEA$MAF <= 2/498/2)) / (nrow(marginEA)+nrow(mafin1000GEA)) # 0.3260734

tmp <- marginEA$BETA; tmp[marginEA$A1 != marginAA$A1] <- - tmp[marginEA$A1 != marginAA$A1]
mean( sign(tmp) == sign(marginAA$BETA) ) # 0.9601227
cor(tmp,marginAA$BETA) # 0.9561637

df.AA <- data.frame(Beta_EA=tmp,Beta_AA=marginAA$BETA,ID=marginAA$ID,stringsAsFactors = F)

df.AA[(df.AA$Beta_AA < -0.5) & (df.AA$Beta_EA > 0.15),] # 893
df.AA[(df.AA$Beta_AA > 0.5) & (df.AA$Beta_EA < -0.2),] # 59
df.AA[(df.AA$Beta_AA > 1) & (df.AA$Beta_EA < 1),] # 710
df.AA[(df.AA$Beta_AA > 0) & (df.AA$Beta_AA < 0.5) & (df.AA$Beta_EA > 1),] # 168

df.AA$eth <- "AA"

df <- rbind(df.EA, df.AA)

write_tsv(df, "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig1.txt")


