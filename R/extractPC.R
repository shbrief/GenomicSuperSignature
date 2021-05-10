#' PCA on gene expression profile
#'
#' Performs a principal components analysis on the given data matrix and returns
#' the results as an object of class \code{\link[stats]{prcomp}}.
#'
#' @param x a numeric or complex matrix (or data frame) which provides the gene
#' expression data for the principal components analysis. Genes in the rows and
#' samples in the columns.
#' @return A result from principle components analysis (PCA) in \code{\link[stats]{prcomp}} object.
#'
#' @examples
#' m = matrix(rnorm(100),ncol=5)
#' extractPC(m)
#'
#' @seealso \code{\link[stats]{prcomp}}
#' @export
extractPC <- function(x) {
    dat.pca <- stats::prcomp(t(x))
    return(dat.pca)
}
