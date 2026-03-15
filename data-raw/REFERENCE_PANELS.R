pacman::p_load(tidyverse, janitor, usethis)
rm(list = ls())

# Build the lookup table for each source ----------------------------------


## RCAv2 -------------------------------------------------------------------

RCA_panels <- RCAv2:::ReferencePanel
RCA_panels[["GlobalPanel_CellTypes"]] <- RCA_panels$GlobalPanel[[1]]
RCA_panels[["GlobalPanel_Tissues"]]   <- RCA_panels$GlobalPanel[[2]]
RCA_panels$GlobalPanel <- NULL
RCA_panels$at <- NULL
RCA_panels$sp <- NULL

levs <- c("GlobalPanel_CellTypes", "GlobalPanel_Tissues",
          "CITESeqPanel", "ENCODEHumanPanel", "ENCODEMousePanel",
          "MonacoPanel", "MonacoBCellPanel", "MonacoMonoPanel", "MonacoTCellPanel",
          "NovershternPanel", "NovershternTCellPanel", "ZhangMouseBrainPanel")

RCAv2_summary <- lapply(RCA_panels, colnames) %>%
  stack() %>%
  dplyr::rename(cell_type = values, panel = ind) %>%
  mutate(tool = "RCAv2", panel = factor(panel, levels = levs)) %>%
  arrange(tool, panel) %>%
  mutate(panel = as.character(panel)) %>%
  select(tool, panel, cell_type)

rm(RCA_panels)

RCAv2_summary %>% pull(panel) %>% unique() %>% length()


## DISCO -------------------------------------------------------------------

DISCO_summary <- DISCO_data("data") %>%
  colnames() %>%
  data.frame(tmp = .) %>%
  mutate(tool = "DISCO") %>%
  separate(tmp, c("cell_type", "panel"), sep = "--") %>%
  select(tool, panel, cell_type)



## SingleR (via celldex) -------------------------------------------------
library(celldex)

SingleR_main_summary <- NULL
SingleR_fine_summary <- NULL

for(ref in listReferences()){

  panel <- fetchReference(ref, "2024-02-26") %>%
    colData() %>%
    data.frame()

  SingleR_main_summary <- bind_rows(
    SingleR_main_summary,
    data.frame(
      tool  = "SingleR",
      panel = paste0(ref, ".main"),
      cell_type = panel$label.main %>% unique() %>% sort())
  )

  SingleR_fine_summary <- bind_rows(
    SingleR_fine_summary,
    data.frame(
      tool  = "SingleR",
      panel = paste0(ref, ".fine"),
      cell_type = panel$label.fine %>% unique() %>% sort())
  )

}


# Combine -----------------------------------------------------------------

allpanels_summary <- bind_rows(RCAv2_summary,
                               DISCO_summary,
                               SingleR_main_summary,
                               SingleR_fine_summary)

dim(allpanels_summary)      # 2194 combinations

allpanels_summary %>%
  select(tool, panel) %>%
  unique() %>%
  filter(!grepl(".fine$", panel)) %>%
  tabyl(tool)
#    tool  n   percent
#   DISCO 39 0.6724138
#   RCAv2 12 0.2068966
# SingleR  7 0.1206897


# Save --------------------------------------------------------------------
use_data(allpanels_summary, overwrite = TRUE, compress = "xz")
