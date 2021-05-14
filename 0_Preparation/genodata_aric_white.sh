#!/usr/bin/env bash
#$ -N aricgeno
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr1.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr1

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr2.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr2

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr3.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr3

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr4.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr4

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr5.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr5

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr6.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr6

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr7.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr7

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr8.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr8

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr9.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr9

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr10.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr10

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr11.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr11

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr12.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr12

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr13.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr13

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr14.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr14

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr15.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr15

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr16.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr16

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr17.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr17

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr18.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr18

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr19.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr19

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr20.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr20

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr21.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr21

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--vcf /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr22.dose.vcf.gz \
--extract-if-info "R2 > 0.8" \
--snps-only \
--maf 0.01 --geno 0.1 --hwe 0.000001 \
--rm-dup exclude-all \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr22


/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr1 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr1

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr2 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr2

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr3 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr3

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr4 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr4

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr5 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr5

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr6 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr6

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr7 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr7

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr8 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr8

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr9 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr9

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr10 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr10

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr11 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr11

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr12 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr12

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr13 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr13

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr14 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr14

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr15 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr15

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr16 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr16

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr17 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr17

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr18 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr18

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr19 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr19

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr20 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr20

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr21 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr21

/users/jzhang2/RESEARCH/tools/plink/plink2 \
--bfile /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/noMatched/white/chr22 \
--update-ids /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/matchID_white.txt \
--remove /dcl01/chatterj/data/jzhang2/eQTLGen/bygene/aricsample_remove.txt \
--make-bed \
--out /dcs01/arking/ARIC_static/ARIC_Data/GWAS/HRC/Aric_HRC_imputation/bedfiles/TOPMed/Filtered/Matched/white/chr22










