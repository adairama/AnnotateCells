pacman::p_load(tidyverse, janitor)

# Load pbmc3k dataset
options(timeout=600)
SeuratData::InstallData("pbmc3k")
pbmc <- SeuratData::LoadData("pbmc3k")
pbmc$orig.ident <- NULL
object.size(pbmc) %>% format(units = "Mb")  # 53.9 Mb

# QC and Filter
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-|^Mt-") %>% round(1)

pbmc@meta.data  <- relocate(pbmc@meta.data,
                            "seurat_annotations",
                            .after = "percent.mt")

VlnPlot(pbmc, c("nCount_RNA", "nFeature_RNA", "percent.mt"))

pbmc <- subset(pbmc,
                 nFeature_RNA > 200 &
                 nFeature_RNA < 2500 &
                 percent.mt   < 5)

pbmc # 2638 cells after QC

## Tabulate the author annotations
pbmc@meta.data %>%
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

# Preprocess the data
pbmc <- pbmc %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

pbmc[["RNA"]]$scale.data <- NULL

# Determine the dimensionality of the data
ElbowPlot(pbmc)
ndims <- 10

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:ndims)

DimPlot(pbmc,
        reduction = "umap",
        group.by  = "seurat_annotations",
        label = TRUE) +
  NoLegend()

object.size(pbmc) %>% format(units = "Mb")  # 55.9 Mb

# Save
pbmc.demo <- pbmc
use_data(pbmc.demo, overwrite = TRUE, compress = "xz")

