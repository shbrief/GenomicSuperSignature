#' z-score each row of the data
#'
#' @param x A gene-expression matrix with genes in rows and samples in columns
#' @author \url{https://github.com/shbrief/PLIER/blob/master/R/Allfuncs.R}
#'
#' @export
rowNorm <- function(x){
  s <- apply(x, 1, sd)
  m <- apply(x, 1, mean);
  x <- sweep(x, 1, m)
  x <- sweep(x, 1, s, "/")
  x
}
