# if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
#
# pacman::p_load( Seurat, SingleR, celldex )
#
# pacman::p_load_gh("prabhakarlab/RCAv2")
# pacman::p_load_gh("lmw123/DISCOtoolkit")
# pacman::p_load_gh("satijalab/Azimuth")
# pacman::p_load_gh("adairama/AnnotateCells")



# Setup -------------------------------------------------------------------
pacman::p_load(tidyverse, Seurat, AnnotateCells)


# Load data ---------------------------------------------------------------
data(pbmc.demo)
pbmc.demo@meta.data %>% head()

pred <- AnnotateCells(pbmc.demo, "RCAv2.GlobalPanel_CellTypes")

pred <- AnnotateCells(pbmc.demo, "DISCO.all")




pbmc.demo@meta.data %>%
  count(seurat_annotations, RCAv2.GlobalPanel_CellTypes) %>%
  arrange(seurat_annotations, desc(n))
