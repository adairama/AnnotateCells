#' Demo dataset prepared from the PBMC3k
#'
#' The famous PBMC3k dataset used in the main Seurat tutorial
#' after quality control filtering and after standard preprocessing:
#' log normalization, highly variable genes, scaling, PCA and UMAP.
#' The dimensionality of the data was set to 10 from Elbow plot inspection.
#'
#' @format A Seurat object with 13,714 genes and 2,635 cells. The meta data slot contains:
#' \describe{
#'   \item{nCount_RNA}{Number of UMIs for each cell}
#'   \item{nFeature_RNA}{Number of genes detected for each cell}
#'   \item{percent.mt}{Percentage of reads mapping to mitochondrial genes per cell}
#'   \item{seurat_annotations}{Cell type annotation provided by authors of Seurat }
#'   \item{RNA_snn_res.0.8}{Communities detected using Louvain algorithm at a resolution 0.8}
#' }
#'
#' @examples
#' pbmc.demo@meta.data %>%
#'   tabyl(seurat_annotations, RNA_snn_res.0.8)
#'
#' DimPlot(pbmc.demo,
#'         reduction = "umap",
#'         group.by  = "RNA_snn_res.0.8",
#'         label = TRUE) +
#'   NoLegend()
#'
#' @source pbmc3k data from the SeuratData package with some data processing
"pbmc.demo"
