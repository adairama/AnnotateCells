if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load( Seurat, SingleR, celldex )

pacman::p_load_gh("prabhakarlab/RCAv2")
pacman::p_load_gh("lmw123/DISCOtoolkit")
pacman::p_load_gh("satijalab/Azimuth")
pacman::p_load_gh("adairama/AnnotateCells")


pacman::p_load(tidyverse, Seurat, AnnotateCells)


data(pbmc.demo)
pbmc.demo@meta.data %>% head()

pred <- AnnotateCells(pbmc.demo, "RCAv2.GlobalPanel_CellTypes")
