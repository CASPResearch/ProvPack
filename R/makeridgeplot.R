#' Make ridge plot of distributional data
#'
#' This function creates overlapping KDEs and plots them as an overlapping 'ridge plot'
#' inspired by the iconic "Unknown Pleasures" album cover by Joy Division. Whilst being very visually
#' attractive it also allows large amounts of KDEs to be plotted in one go, allowing for succint representation.
#' This also allows the different KDE's to be sorted by a clustering algorithm, which groups similar spectra
#' together according to their Kolmogorov-Smirnov statistic.
#'
#' For large datasets this can take a while but be patient!
#' @param distributional An object of class distributional (see "load_distributional"), e.g. Zircon ages, rutile ages etc...
#' @param scale This is the degree of overlap for the ridges. A value of 1 means that the largest overlap in the
#' plot touches the base of the ridge above. Default value is 10, adjust to see what works best for individual cases.
#' @param clustered Boolean, defaults to TRUE indicates whether the samples are to be sorted vertically according
#' to the similarity of their spectra. If FALSE, they are grouped according to their "typing"
#' @param from,to Minimum and maximum x value. Optional
#' @param samebandwidth Boolean, default FALSE. If TRUE, all KDEs are generated with the same bandwidth (not recommended)
#' @param normalise Boolean, default TRUE. All KDEs given same area.
#' @param adaptive Boolean, default TRUE. Geological spectra are often multi-modal, and so the KDE algorithms
#' should have an "adaptive" bandwidth which varies across the x-axis according to point density. This is broad
#' when data is sparse but narrow when data is dense. Turn off only if data is strongly unimodal.
#' @return Returns a ggplot object. Automatically plots it, but can be saved to variable as any normal ggplot object.
#' @keywords distributional, ridges, KDEs
#' @examples
#' \dontrun{ridgeplot <- makeridges(zircs)}
#' \dontrun{ridgeplot <- makeridges(zircs,samebandwidth = FALSE,  normalise = FALSE,  adaptive = FALSE)}
#' \dontrun{ridgeplot <- makeridges(zircs,clustered = FALSE, scale = 20)}
#' \dontrun{ridgeplot <-  makeridges(zircs, from = -500, to = 2500)}
#' \dontrun{ridgeplot <- makeridges(zircs, clustered = FALSE)}
#' @export

makeridges <- function(distributional,
                       scale = 10,
                       clustered = T,
                       from = NA,
                       to = NA,
                       samebandwidth = F,
                       normalise = T,
                       adaptive = T) {
  if(class(distributional)!= "distributional"){
    print("Can only make ridgeplots for distributional objects!")
    return(NULL)
  }
  kdesdf = provenance::KDEs(
    distributional,
    from = from,
    to = to,
    normalise = normalise,
    samebandwidth = samebandwidth,
    adaptive = adaptive
  )

  #Extract into a list of objects KDE
  list_kdes = kdesdf$kdes
  #generate a list and subseqently a dataframe of all the y coordinates
  listys = lapply(list_kdes, function(x)
    x$y)
  df = as.data.frame(listys)
  #extract the common x coordinates
  commonx = list_kdes[[1]]$x
  #add common x coordinates to the dataframe
  df$x = commonx
  #This basically reshapes the data so that x is common, and y's are different. It allows plotting of multiple lines simulataneously
  df_melted <- reshape2::melt(df, id = "x")
  #Adds a typing column by checking which type list each sample name belongs to. Complex code, don't be afraid.
  df_melted$typing = as.character(lapply(df_melted$variable, function(x)
    names(distributional$typing)[indx = which(sapply(distributional$typing, function(y)
      is.element(x, y)))]))
  #Orders the dataframe based on type so they plot in groups next to each other. Manually adjust for visual affect, e.g. by manually setting levels of "sandtype"
  df_melted$typing <-
    factor(df_melted$typing, levels = distributional$TrueOrder)

  if (clustered) {
    ks_dist <- as.dist(provenance::diss(distributional,method='KS'))
    clust <- hclust(ks_dist, method= "ward.D")
    names <- colnames(df)
    clustorder <- names[rev(clust$order)]
    clst <- factor(clustorder, levels = clustorder)
    df_melted$variable <- factor(df_melted$variable, levels = clst)
    sorted <- df_melted
  } else {
    sorted <- df_melted[order(df_melted$typing),]
    sorted$variable <-
      factor(sorted$variable, levels = unique(sorted$variable))
  }
  colnames(sorted) <- c("Variable", "Samples", "Density", "Type")

  plot <- ggplot(data = sorted,
                 aes(
                   x = Variable,
                   y = Samples,
                   height = Density,
                   fill = Type,
                   group = Samples
                 )) +
    ggridges::geom_density_ridges(stat = "identity",
                        scale = scale,
                        size = 0.3) +
    scale_fill_manual(values = as.character(distributional$palette)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0.00, 0)) +
    theme_minimal()
  return(plot)
}

