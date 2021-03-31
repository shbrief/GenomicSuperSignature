#' z-score each row of the data
#'
#' @param x A gene-expression matrix with genes in rows and samples in columns
#' @reference \url{https://github.com/shbrief/PLIER/blob/master/R/Allfuncs.R}
#'
#' @return a matrix with each row mean centered and scaled by rowwise sd
#'
#' @examples
#' x = matrix(rnorm(100),nc=10)
#' y = rowNorm(x)
#' apply(y,1,mean)
#'
#' @export
rowNorm <- function(x){
  return(t(scale(t(x))))
}
