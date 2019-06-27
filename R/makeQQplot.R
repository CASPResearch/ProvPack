#' Make a Lower Quartile - Upper Quartile plot of distributional data
#'
#' This function plots the lower quartile vs the upper quartile for a range of distributions after
#' Anderson et al (2017). Note this is not a typical
#'
#' Andersen, T., Kristoffersen, M., Elburg, M.A., 2017. Visualizing, interpreting and comparing detrital zircon age and Hf isotope data in basin analysis - a graphical approach. Basin Res. https://doi.org/10.1111/bre.12245
#' @param dist An object of class distributional (see "load_distributional"), e.g. Zircon ages, rutile ages etc...
#' @return Returns a ggplot object. Automatically plots it, but can be saved to variable as any normal ggplot object.
#' @keywords distributional, quartile
#' @examples
#' \dontrun{QQplot <- makeQQplot(zircs)}
#' @export
makeQQplot <- function(dist,labels=FALSE) {
  ages = dist$x
  ages_quants <- lapply(ages, quantile)
  names <- names(ages_quants)
  lows <- as.numeric(lapply(ages_quants, function(x)
    x[[2]]))
  upps <- as.numeric(lapply(ages_quants, function(x)
    x[[4]]))

  quant_df = data.frame(names, lows, upps)
  quant_df$typing <-
    as.character(lapply(quant_df$names, function(x)
      names(dist$typing)[indx = which(sapply(dist$typing, function(y)
        is.element(x, y)))]))
  #This sets the order of the sandtypes
  quant_df$typing <-
    factor(quant_df$typing, levels = dist$TrueOrder)

  qq_plot <- ggplot2::ggplot(quant_df,
                    aes(
                      x = lows,
                      y = upps,
                      label = quant_df$names,
                      colour = typing
                    )) +
    geom_point(size = 3) + scale_colour_manual(values = as.character(dist$palette)) +
    theme_minimal() + expand_limits(x = 0, y = 0) +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = "grey",
      size = 1,
      linetype = "dashed"
    ) +
    coord_equal(ratio = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept =
                                                                       0) +
    labs(x = "Lower Quartile", y = "Upper Quartile")
  if(labels){
    qq_plot <- qq_plot + ggrepel::geom_text_repel(size=2.5,show.legend = FALSE)
  }
  return(qq_plot)

}
