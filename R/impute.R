
impute <- function(row, scal, minval) {
  c <- sum(row)
  zeros <- (row == 0)
  replacements <- rep(minval, length(which(zeros))) * scal
  row[zeros] <- replacements
  sumrep <- sum(replacements)

  row[!zeros] <- sapply(row[!zeros], function(x) {
    x <- x * (1 - (sumrep / c))
  })
  return(row)
}
