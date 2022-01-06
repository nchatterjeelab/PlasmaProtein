
## Suppl Fig 2

#### A

library(readr)
library(xlsx)

a <- read_csv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/ExtendedDataFig2.csv")
a <- as.data.frame(a)
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig2.xlsx",
           sheetName = "dat", row.names = FALSE, append = FALSE)

a <- read_csv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/ExtendedDataFig2_leg.csv")
a <- as.data.frame(a)
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig2.xlsx",
           sheetName = "leg", row.names = FALSE, append = TRUE)
