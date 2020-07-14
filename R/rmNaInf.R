#' Remove missing and Inf values from expression dataset
#'
#' @param x A matrix with the expression values.
#' @return The updated input matrix where NA and Inf values are removed.
#'
#' @export
rmNaInf <- function(x) {
  if(!is.matrix(x)) stop("x must be a matrix (expression values)!")
  x <- x[apply(x, 1, function(row) {
    !any(is.na(row) | row %in% c(Inf, -Inf))}), ]
  return(x)
}
