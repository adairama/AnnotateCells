#' Visualize results from one-vs-rest output of FindAllMarkers
#'
#' @param results Output from FindAllMarkers
#' @param order.by
#' @param n
#' @param xthres
#' @param nrow
#' @param text.size
#'
#' @return
#' @export
#'
#' @examples
plot_top_DE <- function(results,
                        order.by = c("avg_log2FC", "pct.diff"),
                        n = 10, xthres = NULL, nrow = NULL, text.size = 3){

  order.by <- match.arg(order.by)
  cat("Ordering by the top", n, "genes based on", order.by, "\n")

  if( order.by == "pct.diff" ){

    if( "pct.diff" %in% colnames(results) ){
      results <- results %>% mutate(score = pct.diff)
    } else {
      warning("pct.diff not detected. Calculating pct.diff = pct.1 - pct.2")
      results <- results %>% mutate(score = pct.1 - pct.2)
    }
  }

  if( order.by == "avg_log2FC" ){
    results <- results %>% mutate(score = avg_log2FC)
  }

  if(is.null(nrow)){
    nrow <- results %>% distinct(cluster) %>% nrow() %>% sqrt() %>% ceiling()
  }

  ## Extract the top genes
  top <- results %>%
    group_by(cluster) %>%
    slice_max(score, n = n) %>%
    arrange(cluster, desc(score)) %>%
    mutate(rank = rank(-score))

  ## Visualization
  if(is.null(xthres)) xthres <- min(top$score)

  g <- ggplot(top, aes(y = rank, x = score, label = gene)) +
    geom_text(size = text.size) +
    facet_wrap( ~ cluster, scales = "free_x", nrow = nrow) +
    labs(x = order.by, y = "Rank", title = NULL) +
    geom_vline(xintercept = xthres, col="red", lty = 2) +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    scale_y_reverse(breaks = unique(top$rank)) +
    theme_bw() +
    theme(panel.grid = element_blank())

  return(g)
}
