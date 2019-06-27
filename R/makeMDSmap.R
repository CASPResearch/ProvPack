#' Make a Multi-Dimensional Scaling (MDS) map of data
#'
#' This is effectively a wrapper for the MDS function in "provenance". Can perform both
#' non-metric and classical MDS. Returns the "stress" of the MDS fit, to assess the goodness of fit.
#'
#' @param obj An object either of class compositional or distributional (see "load_compositional" and "load_distributional")
#' @param classical Boolean, defaults to False. Indicates whether MDS should be non-metric or classical.
#' @param labels Indicates whether labels should be plotted for samples or not. Can be slow for large data-sets...
#' @return Returns a ggplot object. Automatically plots it, but can be saved to variable as any normal ggplot object.
#' @keywords compositional, MDS, distributional
#' @examples
#'\dontrun{makeMDSmap(zircs,labels=T)}
#' @export
makeMDSmap <- function(obj, classical = FALSE,labels=FALSE) {
  if (!(classical)) {
    objMDS <- provenance::MDS(obj)

    #Extracts stress level and saves it as a character string for plotting
    fit <- as.character(round(objMDS$stress, digits = 1))
    fit <- paste("Stress=", fit, sep = " ")
    fit <- paste0(fit, "%")
  }
  else {
    if (obj$method == "SH" | obj$method == "bray") {
      print("Cannot perform classical MDS using Sircombe-Hazleton or Bray dissimilarity.")
      return(NULL)
    } else {
      objMDS <- suppressWarnings(cmdscale(provenance::diss(obj), eig = TRUE))
      eigs <- objMDS$eig
      eigs <- as.vector(sapply(eigs, function(x)
        if (x < 0) {
          x <- 0
        } else {
          x
        }))
      propvar <- cumsum(eigs / sum(eigs))
      fit <-
        paste0(as.character(round(propvar[[2]] * 100, digits = 1)), "%")
      fit <- paste0("Variation explained = ", fit)
    }
  }
  #Creates a "dataframe" of MDS points, with X1=x,X2=y,and X3 = typing
  MDSdf <- data.frame(objMDS$points)
  #This creates a third column to the dataframe with the sandtype of each sample (if available)
  MDSdf$typing <- as.character(lapply(rownames(MDSdf), function(x)
    names(obj$typing)[indx = which(sapply(obj$typing, function(y)
      is.element(x, y)))]))
  MDSdf$typing <- factor(MDSdf$typing, levels = obj$TrueOrder)

  MDSdf$"Sample Name" <- rownames(MDSdf) #Makes a new column in the dataframe so that ggplot can plot them as datapoint labels
  MDSplot <- ggplot(MDSdf,aes(x = X1,y = X2,color = typing,label = MDSdf$"Sample Name")) +
    geom_point(size = 3) + scale_colour_manual(values = as.character(obj$palette)) +
    theme_minimal() + coord_equal(ratio = 1) +
    labs(subtitle = fit,
         x = "Dimension 1",
         y = "Dimension 2")
  if(labels){
    MDSplot <- MDSplot + ggrepel::geom_text_repel(size=2.5)
  }
  return(MDSplot)
}
