#' Remove rows with missing and Inf values from a matrix
#'
#' @param x A numeric matrix.
#' @return The updated input matrix where rows with
#'    NA and Inf values are removed.
#'
#' @examples
#' m = matrix(rnorm(100),ncol=10)
#' m[1,1] = NA
#'
#' m1 = rmNaInf(m)
#' dim(m1)
#'
#' @export
rmNaInf <- function(x) {
  if(!is.matrix(x)) stop("x must be a matrix (expression values)!")
  x <- x[apply(x, 1, function(row) {
    !any(is.na(row) | row %in% c(Inf, -Inf))}), ]
  return(x)
}
