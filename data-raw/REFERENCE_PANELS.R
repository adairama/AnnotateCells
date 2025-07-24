pacman::p_load(tidyverse, janitor, usethis)


# Build the lookup table for each source ----------------------------------

## RCAv2
RCA_panels <- RCAv2:::ReferencePanel
RCA_panels[["GlobalPanel_CellTypes"]] <- RCA_panels$GlobalPanel[[1]]
RCA_panels[["GlobalPanel_Tissues"]]   <- RCA_panels$GlobalPanel[[2]]
RCA_panels$GlobalPanel <- NULL
RCA_panels$at <- NULL
RCA_panels$sp <- NULL

RCAv2_summary <- sapply(RCA_panels, colnames) %>%
  stack() %>%
  rename(cell_type = values, atlas = ind) %>%
  mutate(source = "RCAv2")

rm(RCA_panels)


## DISCO
DISCO_summary <- DISCO_data("data") %>%
  colnames() %>%
  data.frame(tmp = .) %>%
  mutate(source = "DISCO") %>%
  separate(tmp, c("cell_type", "atlas"), sep = "--")

# Consolidate
RCAv2_summary %>% distinct(atlas) %>% nrow()  # 12 atlas
DISCO_summary %>% distinct(atlas) %>% nrow()  # 39 atlas

allpanels_summary <- bind_rows(RCAv2_summary, DISCO_summary)
dim(allpanels_summary)                        # 1501 cell type panels

use_data(allpanels_summary, overwrite = TRUE, compress = "xz")
