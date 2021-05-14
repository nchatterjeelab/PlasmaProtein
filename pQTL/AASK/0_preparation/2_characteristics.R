
library(haven)
#library(data.table)
library(readr)

cov <- read_tsv( paste0('/dcl01/chatterj/data/jzhang2/pwas/pipeline/AASK/aask.cov'))

sum(cov$age_r<=55); mean(cov$age_r<=55)
sum((cov$age_r<=60)&(cov$age_r>=56)); mean((cov$age_r<=60)&(cov$age_r>=56))
sum((cov$age_r<=65)&(cov$age_r>=61)); mean((cov$age_r<=65)&(cov$age_r>=61))
sum(cov$age_r>65); mean(cov$age_r>65)
sum(cov$gender == 2); mean(cov$gender == 2)

