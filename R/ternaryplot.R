#' Makes a ternary plot
#' @param comp An object of type compositional. This should not contain any zeros.
#' @param component1,component2,component3 Strings of the three components to plot against eachother
#' @param centered Boolean indicating whether the ternary phase diagram should be "centered" or not. This
#' is a technique which moves the geometric mean of a data-set to the centre of the plot, thus allowing
#' making variations at the edge of the space visible. see Von Eynatten et al. (2002).
#' @param size Value which adjusts sizes of the points
#' @param isoportions Boolean indicating whether isoPortion lines should be added. Recommended if the data
#' is centered as it allows for the true values to be read off.
#' @return Returns a plot object.
#' @keywords compositional
#' @seealso load_compositional, replacezeros
#' @examples
#' \dontrun{maketernaryplot(qem_rep,"Quartz","Plag","Kfsp",centered=TRUE,size=1,isoportions=TRUE)}
#' @export

maketernaryplot <-
  function(comp,
           component1,
           component2,
           component3,
           centered = FALSE,
           size = 1,
           isoportions = TRUE) {
    vals <- data.frame(comp$x)
    if (any(vals == 0)) {
      print("Contains zeros! Amalgamate variables or replace zeros first and try again")
    } else{
      valscomp <-
        suppressWarnings(compositions::acomp(subset(
          vals, select = c(component1, component2, component3)
        )))
      vals$typing <- as.character(lapply(rownames(vals), function(x)
        names(comp$typing)[indx = which(sapply(comp$typing, function(y)
          is.element(x, y)))]))
      vals$colour <-
        sapply(vals$typing, function(x)
          comp$palette[which(comp$TrueOrder == x)])
      if (isoportions) {
        suppressWarnings(plot.acomp(
          valscomp,
          col = as.character(vals$colour),
          pch = 20,
          center = centered,
          cex = size
        ))
        suppressWarnings(isoPortionLines())
      }
      else{
        suppressWarnings(plot.acomp(
          valscomp,
          col = as.character(vals$colour),
          pch = 20,
          center = centered,
          cex = size
        ))
      }

    }
  }

