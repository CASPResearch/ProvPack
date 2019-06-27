
get_multitypes <- function(df, slist) {
  alltypes <- list()
  for (x in slist) {
    alltypes <- c(alltypes, x$typing)
  }
  out <-
    as.character(lapply(df$SampleNo, function(x)
      names(alltypes)[indx = which(sapply(alltypes, function(y)
        is.element(x, y)))[1]]))
  return(out)
}
#' Makes an INDSCAL map for multiple analyses on the same data-set.
#'
#' See indscal in "provenanace" package documentation, this is effectively just a wrapper.
#' The dissimilarity measures used are inherited from the individual objects.
#'
#' Vermeesch, P., Resentini, A., Garzanti, E., 2016. An R package for statistical provenance analysis. Sediment. Geol., Sediment generation and provenance: processes and pathways 336, 14-25. https://doi.org/10.1016/j.sedgeo.2016.01.009
#' @param ... Multiple objects either of class compositional or distributional or both. For this function to work
#' the data-sets must contain overlapping samples. For example a data-set which has mineralogy, major elements and UPb ages for each sample.
#' @param typing An object of class distributional (see "load_distributional"), e.g. Zircon ages, rutile ages etc...
#' @return Returns a ggplot object. Automatically plots it, but can be saved to variable as any normal ggplot object.
#' @keywords distributional, compositional, INDSCAL
#' @examples
#' \dontrun{inds <- makeindscalmap(zircs,rus,rutmps,typing="typing_palette.csv")}
#' @export
makeindscalmap <- function( ...,typing=NULL,labels=FALSE) {
  ind <- provenance::indscal(...)
  slist <- list(...)
  ind_df <- as.data.frame(ind$gspace)
  ind_df$SampleNo <- rownames(ind_df)
  ind_df$typing = get_multitypes(ind_df, slist)
  if(is.null(typing)){
    palette<-factor(rainbow(length(unique(ind_df$typing))),levels=rainbow(length(unique(ind_df$typing))))
  }else{
  typedf <- read.csv(typing)
  palette <-
    factor(typedf$Colour[which(typedf$Category %in% ind_df$typing)], levels =
             typedf$Colour)
  ind_df$typing <-
    factor(ind_df$typing, levels = typedf$Category)
  }

  ind_plot <- ggplot(data = ind_df,
                     aes(
                       x = D1,
                       y = D2,
                       colour = typing,
                       label = ind_df$SampleNo
                     )) + geom_point(size = 2.5) + scale_colour_manual(values = as.character(palette)) +
    theme_minimal() +
    coord_equal(ratio = 1) +
    labs(x = "Dimension 1",
         y = "Dimension 2",
         caption = "")
  if(labels){
    ind_plot <- ind_plot + ggrepel::geom_text_repel(size = 2.5, show.legend = F)
  }
  xs <- vector()
  ys <- vector()
  dims <- vector()
  for (i in c(1:length(slist))) {
    xs <- c(xs, ind$cweights[[i]][1, 1])
    ys <- c(ys, ind$cweights[[i]][2, 2])
    dims <- c(dims, slist[[i]]$name)
  }

  comp_df <- data.frame(xvals = xs,
                        yvals = ys,
                        compnames = dims)

  ind_plot <-
    ind_plot + ggrepel::geom_text_repel(
      data = comp_df,
      aes(x = xvals * 0.5, y = yvals * 0.5),
      label = comp_df$compnames,
      size = 4,
      colour = "black"
    ) +  geom_segment(
      data = comp_df,
      aes(
        x = 0,
        y = 0,
        xend = xvals * 0.7,
        yend = yvals * 0.7
      ),
      colour = "black",
      arrow = arrow(length = unit(0.4, "cm")),
      inherit.aes = F
    )
  return(ind_plot)
}
#' Makes an INDSCAL map for multiple tests.
#'
#' See procrustes in "provenanace" package documentation, this is effectively just a wrapper.
#' The dissimilarity measures used are inherited from the individual objects.
#'
#' Vermeesch, P., Resentini, A., Garzanti, E., 2016. An R package for statistical provenance analysis. Sediment. Geol., Sediment generation and provenance: processes and pathways 336, 14-25. https://doi.org/10.1016/j.sedgeo.2016.01.009
#' @param ... Multiple objects either of class compositional or distributional or both. For this function to work
#' the data-sets must contain overlapping samples. For example a data-set which has mineralogy, major elements and UPb ages for each sample.
#' @param typing An object of class distributional (see "load_distributional"), e.g. Zircon ages, rutile ages etc...
#' @return Returns a ggplot object. Automatically plots it, but can be saved to variable as any normal ggplot object.
#' @keywords distributional, compositional, procrustes
#' @examples
#' \dontrun{proc <- makeprocrustesmap(zircs,rus,rutmps,typing="typing_palette.csv")}
#' @export
makeprocrustesmap <- function(..., typing=NULL,labels=FALSE) {
  slist <- list(...)
  #Here you now type in the different methods. Zircons = zircs, Apatite = aps etc...
  procru <- provenance::procrustes(...)
  proc_df <- as.data.frame(procru$points)
  proc_df$SampleNo <- procru$labels
  proc_df$typing <- get_multitypes(proc_df, slist)
  if(is.null(typing)){
    proc_palette <- factor(rainbow(length(unique(proc_df$typing))),levels=rainbow(length(unique(proc_df$typing))))
  }else{
  typedf <- read.csv(typing)
  proc_palette <-
    factor(typedf$Colour[which(typedf$Category %in% proc_df$typing)], levels =
             typedf$Colour)
  proc_df$typing <-
    factor(proc_df$typing, levels = typedf$Category)
  }
  proc_plot <- ggplot(data = proc_df,
                      aes(
                        x = V1,
                        y = V2,
                        colour = typing,
                        label = proc_df$SampleNo
                      )) + geom_point(size = 2.5) + scale_colour_manual(values = as.character(proc_palette)) +
    theme_minimal() +
    coord_equal(ratio = 1) +
    labs(x = "Dimension 1",
         y = "Dimension 2")
    if(labels){
      proc_plot <- proc_plot + ggrepel::geom_text_repel(size = 2.5, show.legend = F)
    }
  return(proc_plot)
}
