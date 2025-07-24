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
                panel.data = panel.data[ common, ] %>% as.matrix() ) )
}


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
  projection <- apply(projection, 1, function(x) (x - mean(x))/sd(x) ) %>% t()

  return(projection)
}

RCAv2_annotate <- function(obj, panel.name, ...){

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


#' Load or download and load the DISCO reference data and deg
#'
#' @param type data or deg
#'
#' @return For `type = "data"`, a matrix of weights, where rows are genes and columns are panels. For `type = "deg"`, a dataframe listing the genes involved in each panel.
#' @export
DISCO_data <- function(type = c("data", "deg")){

  ref_path <- file.path( pacman::p_path("AnnotateCells"), "data" )

  # Load reference data. Download first time.
  ref_file <- paste0(ref_path, "/ref_", type, ".rds")

  if (file.exists(ref_file)) {

    message(paste("Loading reference", type, "from local file..."))
    out <- readRDS(ref_file)

  } else {

    message(paste("Downloading reference", type, "from server..."))

    ## Extend timeout temporarily
    op <- options(timeout = 600)
    on.exit(options(op), add = TRUE)

    ## Download
    tryCatch({
      link <- "https://disco.bii.a-star.edu.sg/disco_v3_api/toolkit/getRef"
      if(type == "deg") link <- paste0(link, "Deg")

      con <- curl::curl(link)
      out <- readRDS(con)
      close(con)

      ## Save a local copy
      saveRDS(out, file = ref_file)
      message("Download complete. File saved to: ", ref_file)

    }, error = function(e) {
      stop("Failed to download or read reference data: ", conditionMessage(e))
    })
  }

  return(out)
}


DISCO_annotate <- function(obj, panel.name = "all", assay = "RNA", ...){

  # Prepare the query dataset
  if( nlevels(Idents(obj)) == 1 ) stop("Check if Idents are set")
  Idents(obj) <- paste0("C", Idents(obj))

  data.ave <- AverageExpression(obj, assays = assay, layer = "data")[[1]] %>%
    as.matrix()

  # Prepare the reference dataset
  ref_data <- DISCO_data("data")
  ref_deg  <- DISCO_data("deg")

  # Run group-level prediction
  if(panel.name == "all"){ atlas <- NULL } else { atlas <- panel.name }

  require(DISCOtoolkit)
  grp_pred <- CELLiDCluster(data.ave, ref_data, ref_deg, atlas, ...) %>%
    rownames_to_column("Ident")
  grp_pred %>% mutate(Ident = gsub("^C", "", Ident)) %>% knitr::kable()

  ## Map from group-level to cell-level prediction
  map <- grp_pred %>% select(Ident, predict_cell_type_1) %>% deframe()

  pred <- map[ Idents(obj) ] %>% data.frame()

  dimnames(pred) <- list(colnames(obj),
                         paste0("DISCO.", panel.name))

  return(pred)
}


#' A wrapper function to call various single-cell annotation tools.
#'
#' @param obj A Seurat object
#' @param tool.panel A combination of "toolname.panelname"
#' @param ... Optional parameters
#'
#' @return A data frame with the predicted cell types
#' @export
#'
#' @examples
#' data(pbmc.demo)
#' out <- AnnotateCells( pbmc.demo, "RCAv2.GlobalPanel_CellTypes" )
#'
AnnotateCells <- function(obj, tool.panel, ...){

  ## Check the inputs
  if( !isClass(obj, "Seurat") )
    stop("'obj' needs to be of Seurat class")

  if( identical(GetAssayData(obj, layer = "count"),
                GetAssayData(obj, layer = "data" ) ) )
    warning("The data layer does not appear to be normalized. Please check.")

  tool.name   <- str_split_i(tool.panel, pattern = "\\.", i = 1)
  panel.name  <- str_split_i(tool.panel, pattern = "\\.", i = 2)

  out <- switch(tool.name,
                RCAv2 = RCAv2_annotate(obj, panel.name),
                DISCO = DISCO_annotate(obj, panel.name),
                stop("Unknown tool.name: ", tool.name)
  )

  return(out)
}
