#' Make a principal component biplot of a compositional dataset
#'
#' This function plots both samples and variables (hence the "bi" in biplot) graphically via Principal Component Analysis (PCA).
#' This is effectively a wrapper for the PCA function in "provenance", but visualises it
#' better, especially for large data-sets. Automatically displays components 1 and 2, but
#' can display further ones. Prints out the proportional variation on each component
#' so that it can be seen how many components are significant.
#'
#' Aitchison, J., Greenacre, M., 2002. Biplots of Compositional Data. J. R. Stat. Soc. Ser. C Appl. Stat. 51, 375-392.
#' @param compo An object of class compositional (see "load_compositional"), e.g. mineralogy, major oxides etc... This cannot
#' contain zeros.
#' @param haxis,vaxis These are strings of form "Comp. n" where n is a number. These
#' are by default "Comp.1" and "Comp.2" and they indicate which components to plot on the
#' horizontal and vertical axes respectively. Note that at present there is a bug where the
#' variation displayed on the vertical axis corresponds to the axis which contains less variation
#' irrelevant of whether it is actually the v axis. To be safe, always plot the higher order component
#' on the vertical axis.
#' @param rayscaler Number which adjusts the length of the rays. default is 1, but vary to best
#' fit the data.
#' @return Returns a ggplot object. Automatically plots it, but can be saved to variable as any normal ggplot object.
#' @keywords compositional, biplot, PCA
#' @examples
#' \dontrun{makePCAbiplot(qem_rep,haxis="Comp.1",vaxis="Comp.2")}
#' @seealso maketernaryplot,load_compositional,replacezeros
#' @export
makePCAbiplot <- function(compo,
                          haxis = "Comp.1",
                          vaxis = "Comp.2",
                          rayscaler = 1,
                          labels = FALSE) {
  pca <- provenance::PCA(compo)
  sd_list <- pca$sdev
  indices <-
    which(names(sd_list) == haxis | names(sd_list) == vaxis)

  norm_var <- sd_list ^ 2 / (sum(sd_list ^ 2))
  prop_var <- norm_var / (sum(norm_var))
  print("Cumulate proportional variation:")
  print(cumsum(prop_var))
  tot_var <-
    paste(paste("Total variation =", as.character(round(sum(
      prop_var[indices]
    ), 3) *
      100), sep = " "),
    "%",
    sep = "")
  x_var <-
    paste(paste(" variation contained =", as.character(round(sum(
      prop_var[haxis]
    ), 3) *
      100), sep = " "),
    "%",
    sep = "")
  y_var <-
    paste(paste(" variation contained =", as.character(round(sum(
      prop_var[vaxis]
    ), 3) *
      100), sep = " "),
    "%",
    sep = "")
  pointdata <-
    data.frame(pca$scores,
               samplenames = rownames(pca$scores))
  pointdata$typing <-
    as.character(lapply(rownames(pointdata), function(x)
      names(compo$typing)[indx = which(sapply(compo$typing, function(y)
        is.element(x, y)))]))

  pointdata$typing <-
    factor(pointdata$typing, levels = compo$TrueOrder)
  loadmat <- as.matrix(pca$loadings)
  load_xs <- as.vector(loadmat[, haxis])
  load_ys <- as.vector(loadmat[, vaxis])

  load_df <-
    data.frame(xs = load_xs,
               ys = load_ys,
               vars = rownames(loadmat))
  pcaplot <-
    ggplot(
      data = pointdata,
      aes(
        x = pointdata[, haxis],
        y = pointdata[, vaxis],
        color = typing,
        label = pointdata$samplenames
      )
    ) +
    geom_point(size = 3) + scale_colour_manual(values = as.character(compo$palette)) +
    geom_hline(aes(yintercept = 0), size = 0.2) +
    geom_vline(aes(xintercept = 0), size = 0.2)

  if (labels) {
    pcaplot <- pcaplot + ggrepel::geom_text_repel(size = 2.5)
  }
  pcaplot <- pcaplot +
    suppressWarnings(geom_segment(
      data = load_df,
      aes(
        label = NULL,
        x = 0,
        y = 0,
        xend = xs * 10 * rayscaler,
        yend = ys * 10  * rayscaler
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      alpha = 0.75,
      color = "red"
    )) + ggrepel::geom_text_repel(
      data = load_df,
      aes(
        x = xs * 10 * rayscaler,
        y = ys * 10 * rayscaler,
        label = load_df$vars
      ),
      size = 3.5,
      color = "dark red"
    ) + theme_minimal() + labs(title = tot_var) + xlab(paste(haxis, x_var, sep =
                                                               ",")) + ylab(paste(vaxis, y_var, sep = ","))
  return(pcaplot)
}


