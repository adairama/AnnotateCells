cat0 <- function(...) cat(..., sep = "")

#' Prepare the query and reference data for annotation
#'
#' @description Extracts the data assay from a Seurat object and
#' optionally filters out genes that are expressed in a minimum number of samples.
#' It subsets the filtered data and the reference panel with
#' genes present in both datasets.
#'
#' @param obj Seurat object for annotation
#' @param panel.data Reference panel data
#' @param min.express Filter for genes expressed in at least this
#' proportions of samples in the obj. Expressed between 0 to 1.
#' @param verbose
#'
#' @return A list containing the filtered query and filtered reference panel.
#' @export
data_for_annotation <- function(obj, panel.data, min.express = 0, verbose = TRUE){

  ## Extract the normalized data from Seurat object
  query <- GetAssayData(obj, layer = "data")
  if(verbose) cat0( nrow(query), " genes in query dataset.\n" )

  ## Filter for genes expressed in at least min.express percent of samples
  query <- query[ which(Matrix::rowMeans(query > 0) >= min.express), ]
  if(verbose) cat0( nrow(query), " genes detected in at least ", 100*min.express, "% of the samples.\n" )

  ## Filter for genes overlapping in the reference panel
  if(verbose) cat0(nrow(panel.data), " genes in the reference panel.\n")

  common <- intersect( rownames(query), rownames(panel.data) )
  if(verbose) cat0(length(common), " genes in common used for projection.\n")

  if( length(common) < 30 ) stop("Very few genes detected. Please check the input")

  return( list( query.data = query[ common, ],
                panel.data = panel.data[ common, ] ) )
}


#' Run RCAv2 annotation
#'
#' @param obj Seurat object for annotation
#' @param panel.name Name of the reference panel.
#' @param ... Optional parameters. Currently only `min.express` and `power` is recognized.
#'
#' @return A data frame with the predicted cell types
#' @export
RCAannotate <- function(obj, panel.name, ...){

  ## Get the parameters
  dots <- list(...)
  if( is.null(dots$power) )       dots$power    <- 4
  if( is.null(dots$min.express) ) dots$min.express <- 0.01
  if(dots$min.express < 0 | dots$min.express > 1)
    stop("min.express should be between 0 and 1")

  ## Load and tranform the data
  panel.data <- RCAv2_panels[[panel.name]]
  if(is.null(panel.data)) stop("Check panel name")

  input <- data_for_annotation(obj,
                               panel.data  = panel.data,
                               min.express = dots$min.express )

  ## Find the correlation embedding
  projection <- qlcMatrix::corSparse( input$query.data, input$panel.data )
  dimnames(projection) <- list(colnames(input$query.data),
                               colnames(input$panel.data))

  ## Soft-thresholding and scale each barcode (these are optional actually)
  projection <- abs(projection)^(dots$power) * sign(projection)
  projection <- apply(projection, 1, function(x) (x - mean(x))/sd(x) ) %>% t()

  ## Assign the panel with the highest score as cell type
  pred <- apply(projection, 1, which.max)
  pred <- colnames(projection)[pred]
  pred <- data.frame(pred)
  dimnames(pred) <- list(rownames(projection), paste0("ann.RCAv2.", panel.name))

  return(pred)
}

#' A wrapper function to call various single-cell annotation tools.
#'
#' @param obj A Seurat object
#' @param panel A combination of "toolname.panelname"
#' @param ... Optional parameters
#'
#' @return A data frame with the predicted cell types
#' @export
#'
#' @examples
#' data(pbmc.demo)
#' out <- AnnotateCell( pbmc.demo, "RCAv2.GlobalPanel_CellTypes" )
#'
AnnotateCell <- function(obj,
                         panel = "RCAv2.GlobalPanel_CellTypes",
                         ...){

  if( identical(GetAssayData(obj, layer = "count"), GetAssayData(obj, layer = "data" ) ) )
    warning("The data layer does not appear to be normalized. Please check.")

  tool.name   <- str_split_i(panel, pattern = "\\.", i = 1)
  panel.name  <- str_split_i(panel, pattern = "\\.", i = 2)

  if(tool.name == "RCAv2") out <- RCAannotate(obj, panel.name, ...)

  return(out)
}
