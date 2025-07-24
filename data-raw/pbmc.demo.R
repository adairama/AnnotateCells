pacman::p_load(Seurat, tidyverse, janitor, usethis)


# Load pbmc3k dataset -----------------------------------------------------

options(timeout=600)
SeuratData::InstallData("pbmc3k")
pbmc.demo <- SeuratData::LoadData("pbmc3k")

pbmc.demo$orig.ident <- NULL

pbmc.demo   # 2700 cells
object.size(pbmc.demo) %>% format(units = "Mb")  # 53.9 Mb


# QC and Filter -----------------------------------------------------------

pbmc.demo$percent.mt <- PercentageFeatureSet(pbmc.demo, pattern = "^MT-|^Mt-")

pbmc.demo@meta.data  <- relocate(pbmc.demo@meta.data,
                                 "seurat_annotations",
                                 .after = "percent.mt")

VlnPlot(pbmc.demo, c("nCount_RNA", "nFeature_RNA", "percent.mt"))

pbmc.demo <- subset(pbmc.demo,
                      nFeature_RNA > 200 &
                      nFeature_RNA < 2500 &
                      percent.mt   < 5)

pbmc.demo # 2638 cells after QC



# Tabulate the author annotations -----------------------------------------

pbmc.demo@meta.data %>%
  tabyl(seurat_annotations) %>%
  knitr::kable()

#   |seurat_annotations |   n|   percent|
#   |:------------------|---:|---------:|
#   |Naive CD4 T        | 697| 0.2642153|
#   |Memory CD4 T       | 483| 0.1830933|
#   |CD14+ Mono         | 480| 0.1819560|
#   |B                  | 344| 0.1304018|
#   |CD8 T              | 271| 0.1027293|
#   |FCGR3A+ Mono       | 162| 0.0614102|
#   |NK                 | 155| 0.0587566|
#   |DC                 |  32| 0.0121304|
#   |Platelet           |  14| 0.0053071|



# Preprocess the data -----------------------------------------------------

pbmc.demo <- pbmc.demo %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

pbmc.demo[["RNA"]]$scale.data <- NULL


# Determine the dimensionality of the data
ElbowPlot(pbmc.demo)
ndims <- 10


# Run UMAP
pbmc.demo <- RunUMAP(pbmc.demo, dims = 1:ndims)



# Louvain clustering ------------------------------------------------------

pbmc.demo <- FindNeighbors(pbmc.demo, dims = 1:ndims)
pbmc.demo <- FindClusters (pbmc.demo, res = 0.8)
pbmc.demo$seurat_clusters <- NULL

pbmc.demo$RNA_snn_res.0.8 <-
  formatC( as.numeric(pbmc.demo$RNA_snn_res.0.8),
         format = "d", flag = "0", digits = 1 ) %>%
  paste0("C", .)


# Crosstab with ground truth ----------------------------------------------
pbmc.demo@meta.data %>%
  tabyl(seurat_annotations, RNA_snn_res.0.8) %>%
  knitr::kable()

g1 <- DimPlot(pbmc.demo,
              reduction = "umap",
              group.by  = "seurat_annotations",
              label = TRUE) +
  NoLegend()

g2 <- DimPlot(pbmc.demo,
              reduction = "umap",
              group.by  = "RNA_snn_res.0.8",
              label = TRUE) +
  NoLegend()

g1 | g2


# Save --------------------------------------------------------------------

object.size(pbmc.demo) %>% format(units = "Mb")  # 59.5 Mb

use_data(pbmc.demo, overwrite = TRUE)
