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
pbmc.demo  # 12519 genex x 2638 cells
pbmc.demo@meta.data %>% head()



# Run prediction on selected panels ---------------------------------------
chosen_panels <- c("RCAv2.GlobalPanel_CellTypes",
                   "DISCO.all",
                   "SingleR.hpca.fine")

preds <- AnnotateCells(pbmc.demo, chosen_panels)

pbmc.demo <- AddMetaData(pbmc.demo, preds)


# Visualize the cell level prediction to cluster --------------------------

lapply(to_test, function(combo) {
  align_prediction_to_cluster(
    prediction = norm@meta.data[, combo],
    cluster    = norm$cluster,
    text.size  = 3
  ) + labs(title = combo)
})

pbmc.demo@meta.data %>%
  count(seurat_annotations, RCAv2.GlobalPanel_CellTypes) %>%
  arrange(seurat_annotations, desc(n))
