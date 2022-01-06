
## Suppl Fig 4

#### A

library(xlsx)
library(readr)

a <- as.data.frame(read_csv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/ExtendedDataFig4a.csv")))
write.xlsx(a, file = "//Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_dataExtendedDataFig4.xlsx",
           sheetName = "4a", row.names = FALSE, append = FALSE)

a <- as.data.frame(read_csv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/ExtendedDataFig4b.csv")))
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig4.xlsx",
           sheetName = "4b", row.names = FALSE, append = TRUE)

a <- as.data.frame(read_csv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/ExtendedDataFig4c.csv")))
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig4.xlsx",
           sheetName = "4c", row.names = FALSE, append = TRUE)

a <- as.data.frame(read_csv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/ExtendedDataFig4d.csv")))
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/ExtendedDataFig4.xlsx",
           sheetName = "4d", row.names = FALSE, append = TRUE)




