#' A function to help align prediction labels with cluster labels
#'
#' @param prediction Prediction labels
#' @param cluster Cluster labels
#' @param thres Consider all the prediction labels until thres% of cells are accounted for.
#' @param type The output type:
#' \itemize{
#'  \item{"plot"}{Visual representation of the top prediction labels.}
#'  \item{"split"}{Tabular view of the top prediction labels split by cluster labels.}
#'  \item{"long"}{Tabular view of the top prediction labels.}
#'  \item{"long.all"}{Tabular view of all prediction labels.}
#' }
#'
#' @return Depends on `type` parameter
#' @export
#'
#' @examples
#'
#'
#'
align_prediction_to_cluster <- function(prediction,
                                        cluster,
                                        thres = 0.70,
                                        type = c("plot", "split", "long", "long.all"),
                                        text.size = 2){
  type <- match.arg(type)

  if(!(thres >= 0 & thres <= 1)) stop("`thres` should be between 0 and 1.")

  cts <- data.frame(prediction, cluster) %>%
    count(cluster, prediction) %>%
    arrange(cluster, desc(n)) %>%
    group_by(cluster) %>%
    mutate(prop    = n/sum(n),
           cumprop = cumsum(prop))

  if (type == "long.all") return(cts)

  top_prediction <- cts %>%
    filter(cumprop <= thres | dplyr::lag(cumprop, default = 0) < thres)

  misc <- anti_join(cts,
                    top_prediction %>% select(cluster, prediction),
                    by = c("cluster", "prediction")) %>%
    summarize(n = sum(n), prop = sum(prop)) %>%
    mutate(prediction = "Misc.", cumprop = 1)

  out <- bind_rows(top_prediction, misc) %>%
    arrange(cluster, cumprop)

  rm(cts, top_prediction, misc)

  ## Return outputs
  if (type == "long") return(out)

  if (type == "split") return(split(out %>% as.data.frame(), out$cluster))

  if (type == "plot") {

    tb <- tapply(out$n, out$cluster, sum)

    out$cluster2 <- plyr::mapvalues(out$cluster,
                                    from = names(tb),
                                    paste0(names(tb), " (n=", tb, ")"))

    g <- ggplot(out, aes(x = prop, y = cluster2, fill = prediction)) +
      geom_vline(xintercept = thres) +
      geom_col(position = position_stack(reverse = TRUE)) +
      scale_y_discrete(limits = rev) +
      geom_text(
        aes(label = ifelse(prediction == "Misc.", NA, prediction)),
        position = position_stack(vjust = 0.5, reverse = TRUE),
        size = text.size
      ) +
      scale_x_continuous(
        labels = scales::percent_format(),
        expand = expansion(mult = c(0, 0.01))
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 12, face = "bold")
      )

    cols <- ggplot_build(g)$data[[3]] %>%
      distinct(label, fill) %>%
      deframe()

    cols["Misc."] <- "#D3D3D3"

    g <- g + scale_fill_manual(values = cols)

    return(g)
  }
}
