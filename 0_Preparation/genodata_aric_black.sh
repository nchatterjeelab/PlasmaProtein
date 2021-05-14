#!/usr/bin/env bash
#$ -N aricgeno
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=1000G
#$ -pe local 5
#$ -m e
#$ -M jzhan218@jhu.edu

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr1.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr1

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr2.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr2

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr3.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr3

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr4.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr4

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr5.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr5

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr6.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr6

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr7.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr7

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr8.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr8

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr9.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr9

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr10.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr10

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr11.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr11

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr12.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr12

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr13.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr13

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr14.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr14

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr15.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr15

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr16.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr16

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr17.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr17

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr18.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr18

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr19.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr19

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr20.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr20

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr21.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr21

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/AA/chr22.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr22


/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr1 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr1

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr2 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr2

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr3 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr3

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr4 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr4

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr5 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr5

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr6 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr6

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr7 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr7

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr8 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr8

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr9 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr9

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr10 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr10

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr11 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr11

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr12 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr12

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr13 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr13

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr14 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr14

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr15 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr15

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr16 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr16

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr17 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr17

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr18 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr18

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr19 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr19

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr20 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr20

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr21 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr21

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/black/chr22 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_black.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/black/chr22










