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
  if(!between(min.express, 0, 1)) stop("`min.express` must be between 0 and 1.")
  query <- query[ which(Matrix::rowMeans(query > 0) >= min.express), ]
  if(verbose) cat0( nrow(query), " genes detected in at least ", 100*min.express, "% of the samples.\n" )

  ## Filter for genes overlapping in the reference panel
  if(verbose) cat0(nrow(panel.data), " genes in the reference panel.\n")

  common <- intersect( rownames(query), rownames(panel.data) )
  if(verbose) cat0(length(common), " genes in common used for projection.\n")

  if( length(common) < 30 ) stop("Very few genes detected. Please check the input")

  return( list( query.data = query[ common, ],
                panel.data = panel.data[ common, ] |> as.matrix() ) )
}



# RCA ---------------------------------------------------------------------

#' Get the correlation embeddings / projections for RCAv2 annotation
#'
#' @param obj Seurat object for annotation
#' @param panel.name Name of the reference panel.
#' @param ... Optional parameters. Currently only `min.express` and `power` is recognized.
#'
#' @return A data frame with the predicted cell types
#' @export
RCA_project <- function(obj, panel.name = NULL, ...){

  dots <- list(...)
  if( is.null(dots$power) )       dots$power       <- 4
  if( is.null(dots$min.express) ) dots$min.express <- 0.01
  if(dots$min.express < 0 | dots$min.express > 1)
    stop("min.express should be between 0 and 1")

  ## Load the reference panel
  RCA_panels <- RCAv2:::ReferencePanel

  if( panel.name == "GlobalPanel_CellTypes" ) {

    panel.data <- RCA_panels[["GlobalPanel"]][[1]]
    panel.data[ panel.data < 1.5 ] <- 1.5

  } else if( panel.name == "GlobalPanel_Tissues" ) {

    panel.data <- RCA_panels[["GlobalPanel"]][[2]]
    panel.data[ panel.data < 0.75 ] <- 0.75

  } else {

    if( !(panel.name %in% names(RCA_panels)) )
      stop("Invalid panel.name: ", panel.name)

    panel.data <- RCA_panels[[panel.name]]

  }

  ## Subset
  cat("\n", panel.name, "\n")
  input <- data_for_annotation(obj,
                               panel.data  = panel.data,
                               min.express = dots$min.express )

  ## Find the correlation embedding
  projection <- qlcMatrix::corSparse( input$query.data, input$panel.data )
  dimnames(projection) <- list(colnames(input$query.data),
                               colnames(input$panel.data))

  ## Soft-thresholding and scale each barcode
  projection <- abs(projection)^(dots$power) * sign(projection)
  projection <- apply(projection, 1, function(x) (x - mean(x))/sd(x) ) |> t()

  return(projection)
}

RCAv2_annotate <- function(obj, panel.name, ...){

  require(RCAv2)

  if( panel.name == "GlobalPanel" ) {

    proj <- cbind(
      RCA_project(obj, "GlobalPanel_CellTypes", ...),
      RCA_project(obj, "GlobalPanel_Tissues", ...)
    )

  } else {

    proj <- RCA_project(obj, panel.name, ...)

  }

  ## Assign the panel with the highest score as cell type
  pred <- apply(proj, 1, which.max)
  pred <- colnames(proj)[pred]
  pred <- data.frame(pred)
  dimnames(pred) <- list(rownames(proj),
                         paste0("RCAv2.", panel.name))

  return(pred)
}


# DISCO -------------------------------------------------------------------
DISCO_annotate <- function(obj, panel.name = "all", local = FALSE,
                           n_predict = 1, ref_path = "~/.disco_cache", verbose = TRUE){

  # Prepare the query dataset
  if( nlevels(Idents(obj)) == 1 ) stop("Check if Idents are set")

  avg <- Seurat::AverageExpression(obj, assays = "RNA", layer = "counts")$RNA |>
    as.matrix()

  colnames(avg) <- gsub("^g(\\d+)$", "\\1", colnames(avg))

  ## Calculate
  require(DISCOtoolkit)

  if (local) {
    grp_pred <- DISCOtoolkit:::.cellid_local (avg, panel.name, n_predict, ref_path, verbose)
  } else {
    grp_pred <- DISCOtoolkit:::.cellid_remote(avg, panel.name, n_predict, verbose)
  }

  grp_pred |> knitr::kable() |> print()


  ## Map from group-level to cell-level prediction
  map <- grp_pred |> select(cluster, predicted_cell_type) |> deframe()

  pred <- map[ Idents(obj) ] |> data.frame()

  dimnames(pred) <- list(colnames(obj),
                         paste0("DISCO.", panel.name))

  return(pred)
}


# SingleR -----------------------------------------------------------------

#' Returns the pruned labels from SingleR
#'
#' @param obj A Seurat object
#' @param panel.name Reference panel. Run `celldex::listReferences()` to see available panels.
#' @param version Currently fixed to "2024-02-26" and can be updated if there are new releases.
#' @param label Either `label.main`, `label.fine` or `label.ont`. If NULL, it will revert to label.main
#' @param ...
#'
#' @return A data frame that returns the `pruned.labels` column of SingleR
#' @export
#'
#' @examples
SingleR_annotate <- function(obj, panel.name, version = "2024-02-26", label = NA, ...){

  label <- ifelse(is.na(label), "main", label)
  label <- paste0("label.", label)
  stopifnot( label %in% c("label.main", "label.fine", "label.ont") )

  ref.se <- celldex::fetchReference(panel.name, version, realize.assays = TRUE)

  if( suppressWarnings( nrow(GetAssayData(obj, layer = "scale.data")) == 0 ) ){
    obj <- obj |> ScaleData(verbose = FALSE)
  }

  suppressMessages(
    pred <- SingleR::SingleR(test = as.SingleCellExperiment(obj),
                             ref  = ref.se,
                             labels = SummarizedExperiment::colData(ref.se)[ , label])
  )

  ## Return the pruned labels
  pred <- pred |>
    data.frame() |>
    select(y = pruned.labels) |>
    mutate(y = ifelse(is.na(y), "Unassigned", y))

  colnames(pred) <- paste0("SingleR.", panel.name,
                           gsub("label", "", label) )
  return(pred)

}


#' A wrapper function to call one annotation tools.
#'
#' @param obj A Seurat object
#' @param combo combination of "toolname.panelname"
#' @param label Only used for SingleR too. Values are either `label.main`, `label.fine` or `label.ont`. If NULL, it will revert to label.main
#' @param ... Optional parameters
#'
#' @return A data frame with the predicted cell types
#' @export
#'
#' @examples
#' data(pbmc.demo)
#' out <- AnnotateCells( pbmc.demo, "RCAv2.GlobalPanel_CellTypes" )
#'
AnnotateCells_internal <- function(obj, combo, label = NULL, ...){

  cat("\n", combo)

  tool.name  <- str_split_i(combo, pattern = "\\.", i = 1)
  panel.name <- str_split_i(combo, pattern = "\\.", i = 2)
  label.name <- str_split_i(combo, pattern = "\\.", i = 3)

  out <- switch(tool.name,
                RCAv2   = RCAv2_annotate(obj, panel.name, ...),
                DISCO   = DISCO_annotate(obj, panel.name, ...),
                SingleR = SingleR_annotate(obj, panel.name,
                                           label = label.name, ...),
                stop(tool.name, "not available.")
  )

  return(out)
}


#' A wrapper function to call multiple annotation tools.
#'
#' @param obj A Seurat object
#' @param combos Combinations of "toolname.panelname"
#' @param label Only used for SingleR too. Values are either `label.main`, `label.fine` or `label.ont`. If NULL, it will revert to label.main
#' @param ... Optional parameters
#'
#' @return A data frame with the predicted cell types
#' @export
#'
#' @examples
#' data(pbmc.demo)
#' out <- AnnotateCells( pbmc.demo, "RCAv2.GlobalPanel_CellTypes" )
#'
AnnotateCells <- function(obj, combos, label = NULL, ...){

  require(dplyr)
  require(stringr)

  ## Check the inputs
  if( !isClass(obj, "Seurat") )
    stop("'obj' needs to be of Seurat class")

  if( identical(GetAssayData(obj, layer = "count"),
                GetAssayData(obj, layer = "data" ) ) )
    warning("The data layer does not appear to be normalized. Please check.")

  ## Call the internal function iteratively
  preds <- lapply(combos, AnnotateCells_internal, obj = obj)
  preds <- do.call(cbind, preds)

  return(preds)
}
