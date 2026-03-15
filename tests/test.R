# if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
#
# pacman::p_load( Seurat, SingleR, celldex )
# pacman::p_load_gh("prabhakarlab/RCAv2")
# pacman::p_load_gh("lmw123/DISCOtoolkit")
# pacman::p_load_gh("satijalab/Azimuth")
# pacman::p_load_gh("adairama/AnnotateCells")





# Setup -------------------------------------------------------------------
pacman::p_load(tidyverse, Seurat, AnnotateCells)
rm(list = ls())


# Load data ---------------------------------------------------------------
data(pbmc.demo)
pbmc.demo@meta.data %>% head()



# Run prediction on selected panels ---------------------------------------
pred1 <- AnnotateCells(pbmc.demo, "RCAv2.GlobalPanel_CellTypes")
pred2 <- AnnotateCells(pbmc.demo, "DISCO.all")
pred3 <- AnnotateCells(pbmc.demo, "SingleR.hpca")



pbmc.demo@meta.data %>%
  count(seurat_annotations, RCAv2.GlobalPanel_CellTypes) %>%
  arrange(seurat_annotations, desc(n))
