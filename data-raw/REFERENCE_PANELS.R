pacman::p_load(tidyverse, janitor, usethis)

# RCAv2 panels ----
RCAv2_panels <- RCAv2:::ReferencePanel

tmp <- RCAv2_panels$GlobalPanel[[1]]
tmp[ tmp < 1.5 ] <- 1.5
RCAv2_panels[["GlobalPanel_CellTypes"]] <- tmp
rm(tmp)

tmp <- RCAv2_panels$GlobalPanel[[2]]
tmp[ tmp < 0.75 ] <- 0.75
RCAv2_panels[["GlobalPanel_Tissue"]] <- tmp
rm(tmp)

RCAv2_panels$GlobalPanel <- NULL
RCAv2_panels$at <- NULL
RCAv2_panels$sp <- NULL

RCAv2_panels <- lapply(RCAv2_panels, function(x) as.matrix(x))

use_data(RCAv2_panels, overwrite = TRUE, compress = "xz")
rm(RCAv2_panels)


# DISCO panels ----
library(DISCOtoolkit)
options(timeout = 600)

## DISCO DEG ----
DISCO_ref_deg <- readRDS( url(paste0(getOption("disco_url"), "toolkit/getRefDeg")) )
dim(DISCO_ref_deg)  # 197,516 genes

DISCO_ref_deg %>% pull(group) %>% unique() %>% length()  # 1098 cell types

DISCO_ref_deg %>% tabyl(group) %>% arrange(n) %>% pull(n) %>% summary()
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    3.0    84.0   147.0   179.9   232.0  1009.0

use_data(DISCO_ref_deg, overwrite = TRUE, compress = "xz")

## DISCO reference panel ----
DISCO_ref_data <- readRDS( url(paste0(getOption("disco_url"), "toolkit/getRef")) )

dim(DISCO_ref_data)
# 33538 genes x 1098 cell tpes

use_data(DISCO_ref_data, overwrite = TRUE, compress = "xz")


## DISCO atlas ----
DISCO_atlas <- str_split_i( colnames(DISCO_ref_data), pattern = "--", i = 2 ) %>% unique()
length(DISCO_atlas) # 39 atlas

use_data(DISCO_atlas, overwrite = TRUE, compress = "bzip2")

rm(DISCO_ref_deg, DISCO_ref_data, DISCO_atlas)
