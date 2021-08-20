
rm(list = ls())

for (chr in 1:22){
    suppressMessages(system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2",
                                   " --bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/White/update_rsid/chr",chr,
                                   " --keep /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pheno.id",
                                   " --make-bed",
                                   " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/imputation_expression_on_aric/geno/chr",chr)))
}