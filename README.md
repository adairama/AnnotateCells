
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AnnotateCell

<!-- badges: start -->
<!-- badges: end -->

Cell type identification remains a key challenge in single-cell RNA-seq
and spatial transcriptomics datasets despite a decade of research in
this area.

There are many reference-based annotation tools to help with cell type
annotation. However, the input and output formats differ greatly,
increasing the burden on users to use and compare multiple tools.

This package `AnnotateCell` aims to provide a unified wrapper to run
several annotation tools conveniently.

## Installation

You can install AnnotateCell like so:

``` r
library(remotes)
install_github("adairama/AnnotateCell")
```

# Illustration

## Demo dataset

We will load the `pbmc.demo` dataset. This is a Seurat object containing
2,635 peripheral blood mononuclear cells (PBMC) from 10X Genomics
experiment. This is the same dataset used in the fundamental Seurat
vignette after QC filtering with PCA and UMAP embedding. See
help(pbmc.demo) for more details.

``` r
library(AnnotateCell)
#> Loading required package: Seurat
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
#> Loading required package: tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

pbmc.demo
#> An object of class Seurat 
#> 13714 features across 2635 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, umap

pbmc.demo@meta.data %>% head()
#>                nCount_RNA nFeature_RNA percent.mt seurat_annotations
#> AAACATACAACCAC       2419          779        3.0       Memory CD4 T
#> AAACATTGAGCTAC       4903         1352        3.8                  B
#> AAACATTGATCAGC       3147         1129        0.9       Memory CD4 T
#> AAACCGTGCTTCCG       2639          960        1.7         CD14+ Mono
#> AAACCGTGTATGCG        980          521        1.2                 NK
#> AAACGCACTGGTAC       2163          781        1.7       Memory CD4 T
```

The authors of Seurat package provided their annotation which is stored
as `seurat_annotations` in the meta data.

``` r
pbmc.demo@meta.data %>% 
  janitor::tabyl(seurat_annotations)
#>  seurat_annotations   n     percent
#>         Naive CD4 T 697 0.264516129
#>        Memory CD4 T 483 0.183301708
#>          CD14+ Mono 479 0.181783681
#>                   B 344 0.130550285
#>               CD8 T 269 0.102087287
#>        FCGR3A+ Mono 162 0.061480076
#>                  NK 155 0.058823529
#>                  DC  32 0.012144213
#>            Platelet  14 0.005313093
```

## Run RCAv2 annotation

Here is an example on how to run the RCAv2 tool with the
GlobalPanel_CellTypes panel.

``` r
pred <- AnnotateCell(pbmc.demo, "RCAv2.GlobalPanel_CellTypes")
#> 13714 genes in query dataset.
#> 8882 genes detected in at least 1% of the samples.
#> 5209 genes in the reference panel.
#> 1192 genes in common used for projection.

head(pred)
#>                ann.RCAv2.GlobalPanel_CellTypes
#> AAACATACAACCAC    L74_T.Cell_CD4.Centr..Memory
#> AAACATTGAGCTAC          L51_B.Cell_Bone.Marrow
#> AAACATTGATCAGC    L75_T.Cell_CD4.Centr..Memory
#> AAACCGTGCTTCCG               L60_Monocyte_CD16
#> AAACCGTGTATGCG              L86_NK.Cell_CD56Lo
#> AAACGCACTGGTAC    L75_T.Cell_CD4.Centr..Memory
```

## Comparing the annotation

The first task is to add the prediction to the Seurat object.

``` r
pbmc.demo <- AddMetaData(pbmc.demo, pred)
```

``` r
cts <- pbmc.demo@meta.data %>% 
  count(group = seurat_annotations, pred = ann.RCAv2.GlobalPanel_CellTypes) %>% 
  arrange(group, desc(n))
```
