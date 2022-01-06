
## Fig 2

#### A

library(xlsx)

a <- read.table(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/2a.txt"))
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig2.xlsx",
           sheetName = "2a", row.names = FALSE, append = FALSE)

a <- read.table(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/2b.txt"))
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig2.xlsx",
           sheetName = "2b", row.names = FALSE, append = TRUE)

a <- read.table(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/2c.txt"), header = T)
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig2.xlsx",
           sheetName = "2c", row.names = FALSE, append = TRUE)

a <- read.table(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/2d.txt"), header = T)
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig2.xlsx",
           sheetName = "2d", row.names = FALSE, append = TRUE)

a <- read.table(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/2e.txt"), header = T)
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig2.xlsx",
           sheetName = "2e", row.names = FALSE, append = TRUE)

a <- read.table(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/500Kb/*Source_data/preparation/2f.txt"), header = T)
write.xlsx(a, file = "/Users/jnz/Dropbox/PWAS_manuscript/NatureGenetics/2021_12_revision4/Final_files_prepared_for_submission/Source_data/Fig2.xlsx",
           sheetName = "2f", row.names = FALSE, append = TRUE)


