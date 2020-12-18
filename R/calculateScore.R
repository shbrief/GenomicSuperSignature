#' Calculate the validation score for a new dataset
#'
#' @import methods
#'
#' @param dataset A gene expression matrix, where genes are in rows and rownames
#' are in 'symbol' format. It can be SummarizedExperiment, ExpressionSet, or matrix objects.
#' @param RAVmodel PCAGenomicSignatures object. Output from \code{buildAvgLoading}
#' function, a matrix of avaerage loadings, can be directly provided.
#' @param rescale.after If it is \code{TRUE}, the continuous scores are rescaled
#' post assignment, so avgerage loadings have the same standard deviation in different
#' studies. If it is \code{FALSE}, the rescaling of column (= dividing by \code{sqrt(sum(x^2})
#' is done before score assignment. Default is \code{TRUE}.
#'
#' @return A list containing the score matrices for input datasets. Scores are
#' assigned to each sample (row) for each cluster (column).
#'
#' @examples
#' data(miniRAVmodel)
#' library(bcellViper)
#' data(bcellViper)
#' score <- calculateScore(dset, miniRAVmodel)
#'
#' @export
calculateScore <- function(dataset, RAVmodel, rescale.after = TRUE) {

    if (is(RAVmodel, "PCAGenomicSignatures")) {
        avg.loadings <- assay(RAVmodel)
    } else if (class(RAVmodel) %in% c("data.frame", "matrix")) {
        avg.loadings <- RAVmodel
    }

    if (!is.list(dataset)) {dataset <- list(dataset)}
    validationScore <- lapply(dataset, function(dat) {
        if (is(dat, "ExpressionSet")) {
            count <- Biobase::exprs(dat)
        } else if (is(dat, "SummarizedExperiment")) {
            count <- SummarizedExperiment::assay(dat)
        } else if (is(dat, "matrix")) {
            count <- dat
        }

        count <- count[apply(count, 1, function(x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]
        count <- apply(count, 1, function(x) {x - mean(x)}) %>% t  # row centered
        gene_common <- intersect(rownames(avg.loadings), rownames(count))

        if (isFALSE(rescale.after)) {
          score <- t(count[gene_common,]) %*% apply(avg.loadings[gene_common,], 2,
                                                    function(x) x / sqrt(sum(x^2, na.rm = TRUE)))
        } else {
          # CRC paper version
          score <- t(count[gene_common,]) %*% as.matrix(avg.loadings[gene_common,])
          score <- (t(score) / apply(score, 2, stats::sd)) %>% t
        }

        colnames(score) <- colnames(avg.loadings)
        return(score)
    })

    if (length(dataset) == 1) {return(validationScore[[1]])} else {return(validationScore)}
}
