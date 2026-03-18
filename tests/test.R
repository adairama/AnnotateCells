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


pbmc.demo <- AddMetaData(pbmc.demo, pred2)

align_prediction_to_cluster(
  pbmc.demo$RCAv2.GlobalPanel_CellTypes,
  pbmc.demo$RNA_snn_res.0.8,
  text.size = 6
)

align_prediction_to_cluster(
  pbmc.demo$DISCO.all,
  pbmc.demo$RNA_snn_res.0.8,
  text.size = 6
)


align_prediction_to_cluster(
  pbmc.demo$SingleR.hpca.fine,
  pbmc.demo$RNA_snn_res.0.8,
  text.size = 6
)

pbmc.demo@meta.data %>%
  count(seurat_annotations, RCAv2.GlobalPanel_CellTypes) %>%
  arrange(seurat_annotations, desc(n))
