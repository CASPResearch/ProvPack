#' Replaces zeros in a compositional dataset by multiplicative replacement
#'
#' Log ratio analysis of compostiional data makes zero values in data problematic
#' because the logarithm of a zero is undefined. Thus, if these are treated as rounded
#' zeros they should be replaced. However, if more than 10% of the data-set is zeros it
#' this can be problematic. This function replaces the zeros in a data-set with a small
#' value equal to 0.65 times either the detection limit, or, if this is not provided
#' the minimum non-zero value in the data-set. The stategy used is outlined in Martin-Fernandez et al. 2003.
#' It is strongly recommended that sensitivity tests should be performed by varying the value of "scal".
#'
#' Martín-Fernández, J.A., Barceló-Vidal, C., Pawlowsky-Glahn, V., 2003. Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation. Math. Geol. 35, 253–278. https://doi.org/10.1023/A:1023866030544
#' @seealso fractionzeros
#' @param comp Object of type compositional
#' @param scal Scalar value, default at 0.65, but can be varied for sensitivity tests
#' @param detlim The minimum detection limit for data. If not specified, it is approximated as the smallest nonzero value.
#' @keywords compositional, zeros, zero-replacement
#' @export
replacezeros <- function(comp, scal = 0.65, detlim = NULL) {
  vals <- comp$x
  minval <- vector()
  for (i in colnames(vals)) {
    minval <- min(c(minval, unique(vals[order(vals[, i]), i])[[2]]))
  }
  if (is.null(detlim)) {
    detlim <- minval
  }

  for (j in rownames(vals)) {
    vals[j,] <- impute(vals[j,], scal, detlim)
  }
  comp$x <- vals
  return(comp)
}
