.validateScore <- function(dat, rescale.after, avg.loadings) {
    count <- .extractExprsMatrix(dat)
    count <- count[apply(count, 1,
                         function(x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]
    count <- apply(count, 1, function(x) {x - mean(x)}) %>% t  # row centered
    gene_common <- intersect(rownames(avg.loadings), rownames(count))

    if (rescale.after) {
        # CRC paper version
        score <- t(count[gene_common,]) %*% as.matrix(avg.loadings[gene_common,])
        score <- (t(score) / apply(score, 2, stats::sd)) %>% t
    } else {
        score <- t(count[gene_common,]) %*%
            apply(avg.loadings[gene_common,], 2,
                  function(x) x / sqrt(sum(x^2, na.rm = TRUE)))
    }

    colnames(score) <- colnames(avg.loadings)
    return(score)
}

#' Calculate the validation score for a new dataset
#'
#' @import methods
#'
#' @param dataset A gene expression dataset to validate. It can be ExpressionSet,
#' SummarizedExperiment, RangedSummarizedExperiment, or matrix. Rownames (genes) should be in
#' 'symbol' format. If it is a matrix, genes should be in rows and samples in columns.
#' @param RAVmodel PCAGenomicSignatures object. Output from \code{buildAvgLoading}
#' function, a matrix of average loadings, can be directly provided.
#' @param rescale.after If it is \code{TRUE}, the continuous scores are rescaled
#' post assignment, so average loadings have the same standard deviation in different
#' studies. If it is \code{FALSE}, the rescaling of column (= dividing by \code{sqrt(sum(x^2})
#' is done before score assignment. Default is \code{TRUE}.
#'
#' @return A list containing the score matrices for input datasets. Scores are
#' assigned to each sample (row) for each cluster (column).
#'
#' @examples
#' data(miniRAVmodel)
#' dim(miniRAVmodel)
#' miniRAVmodel
#' if(require(bcellViper)) {
#'   data(bcellViper)
#'   dset
#'   head(rownames(dset))
#'   score <- calculateScore(dset, miniRAVmodel)
#'   dim(score)
#'   head(colnames(score))
#'   head(rownames(score))
#' }
#' @export
calculateScore <- function(dataset, RAVmodel, rescale.after = TRUE) {

    if (is(RAVmodel, "PCAGenomicSignatures")) {
        avg.loadings <- assay(RAVmodel)
    } else if (is.data.frame(RAVmodel) | is.matrix(RAVmodel)) {
        avg.loadings <- RAVmodel
    }

    if (!is.list(dataset)) {
        validationScore <- .validateScore(dataset, rescale.after, avg.loadings)
    } else {
        validationScore <- lapply(dataset, .validateScore,
                                  rescale.after, avg.loadings)
    }

    return(validationScore)
}
