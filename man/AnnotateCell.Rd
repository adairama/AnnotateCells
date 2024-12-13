% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation_functions.R
\name{AnnotateCell}
\alias{AnnotateCell}
\title{A wrapper function to call various single-cell annotation tools.}
\usage{
AnnotateCell(obj, panel = "RCAv2.GlobalPanel_CellTypes", ...)
}
\arguments{
\item{obj}{A Seurat object}

\item{panel}{A combination of "toolname.panelname"}

\item{...}{Optional parameters}
}
\value{
A data frame with the predicted cell types
}
\description{
A wrapper function to call various single-cell annotation tools.
}
\examples{
data(pbmc.demo)
out <- AnnotateCell( pbmc.demo, "RCAv2.GlobalPanel_CellTypes" )

}
