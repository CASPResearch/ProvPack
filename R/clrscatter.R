#' Plots a scattergram of two clr transformed variables.
#'
#' The nature of compositional data means that plotting raw percentages against eachother on a scatterplot
#' is fraught with problems. This can be improved by transforming the variables into clr space. Ignore the warnings.
#'
#' Pearson, K., 1897. Mathematical contributions to the theory of evolution.—On a form of spurious correlation which may arise when indices are used in the measurement of organs. Proc. R. Soc. Lond. 60, 489–498. https://doi.org/10.1098/rspl.1896.0076
#' @seealso clr, makePCAbiplot, load_compositional
#' @param comp Object of type compositional
#' @param xval,yval String referring to the names of the variables on x and y axes
#' @examples
#' \dontrun{clrscatterplot(qem_rep,"Diagenetic","Porosity")}
#' @keywords compositional, clr, scatterplot
#' @export
clrscatterplot <- function(comp, xval, yval,labels=FALSE) {
  vals <- suppressWarnings(compositions::acomp(comp$x))
  clrvals <- suppressWarnings(compositions::clr(vals))
  xs <- clrvals[, xval]
  ys <- clrvals[, yval]
  if (any(xs == 0) | any(ys == 0)) {
    print("Warning: zero-values present, replacement recommended")
  }
  names <- rownames(clrvals)
  typing <- as.character(lapply(names, function(x)
    names(comp$typing)[indx = which(sapply(comp$typing, function(y)
      is.element(x, y)))]))
  typing <- factor(typing, levels = comp$TrueOrder)
  df <- data.frame(
    x = xs,
    y = ys,
    SampleNo = names,
    Type = typing
  )
  scatterplot <-
    ggplot(data = df, aes(
      x = x,
      y = y,
      colour = Type,
      label = SampleNo
    )) + geom_point(size = 2) +
    scale_colour_manual(values = as.character(comp$palette)) + theme_minimal()  + coord_equal() +
    xlab(paste0("c", "(", xval, ")")) + ylab(paste0("c", "(", yval, ")"))
  if(labels){
    scatterplot <- scatterplot+ ggrepel::geom_text_repel()
  }
  return(scatterplot)
}

