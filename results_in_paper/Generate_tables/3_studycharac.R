
####################################
####################################
## pQTL summary

## white

rm(list=ls())

library(readr)
library(stringr)

aric.cov <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/aric.cov', header=TRUE, stringsAsFactors=F)
dim(aric.cov)
summary(aric.cov$v3age31)
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  #50.00   55.00   60.00   60.33   65.00   72.00

sum(aric.cov$gender==1) # 3429
sum(aric.cov$gender==2) # 3784
sum(aric.cov$v3center=="M") # 2707
sum(aric.cov$v3center=="F") # 2187
sum(aric.cov$v3center=="W") # 2319

sum(aric.cov$v3age31<=55)
sum((aric.cov$v3age31>55)&(aric.cov$v3age31<=60))
sum((aric.cov$v3age31>60)&(aric.cov$v3age31<=65))
sum(aric.cov$v3age31>65)

aric.cov <- read.table('/dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/aric.cov', header=TRUE, stringsAsFactors=F)
dim(aric.cov)
summary(aric.cov$v3age31)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  #49.00   54.00   59.00   59.18   64.00   73.00
sum(aric.cov$gender==1) # 715
sum(aric.cov$gender==2) # 1156
sum(aric.cov$v3center=="J") # 1645
sum(aric.cov$v3center=="F") # 226

sum(aric.cov$v3age31<=55)
sum((aric.cov$v3age31>55)&(aric.cov$v3age31<=60))
sum((aric.cov$v3age31>60)&(aric.cov$v3age31<=65))
sum(aric.cov$v3age31>65)


