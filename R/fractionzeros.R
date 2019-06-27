#' Returns the fraction of a compositional data set that is zeros
#'
#' Log ratio analysis of compostiional data makes zero values in data problematic
#' because the logarithm of a zero is undefined. Thus, if these are treated as rounded
#' zeros they should be replaced. However, if more than 10% of the data-set is zeros it
#' this can be problematic. This function returns the proportion of the data-set that is zeros
#' so the user can assess whether amalgamation or replacement of zeros should go ahead
#' If there are zeros in a data-set, PCA cannot be performed.
#'
#' Martín-Fernández, J.A., Barceló-Vidal, C., Pawlowsky-Glahn, V., 2003. Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation. Math. Geol. 35, 253–278. https://doi.org/10.1023/A:1023866030544
#' @seealso replacezeros
#' @param comp Object of type compositional
#' @keywords compositional, zeros, zero replacement
#' @export
#'
fractionzeros <- function(comp) {
  df <- comp$x
  num_zeros <- sum(sapply(df, function(x)
    length(which(x == 0))))
  num_vals <- nrow(df) * ncol(df)
  return(num_zeros / num_vals)
}