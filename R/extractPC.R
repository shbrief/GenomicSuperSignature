#' PCA on gene expression profile
#'
#' Performs a principal components analysis on the given data matrix and returns
#' the results as an object of class prcomp.
#'
#' @param x a numeric or complex matrix (or data frame) which provides the gene
#' expression data for the principal components analysis. Genes in the rows and
#' samples in the columns.
#' @return A result from principle components analysis (PCA) in \code{\link[stats]{prcomp}} object.
#'
#' @seealso \code{\link[stats]{prcomp}}
#' @export
extractPC <- function(x) {
    dat.pca <- prcomp(t(x))
    return(dat.pca)
}
