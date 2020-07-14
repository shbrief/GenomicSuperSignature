#' Calculate the validation score for a new dataset
#'
#' @import methods
#'
#' @param dataset A gene expression matrix, where genes are in rows and rownames
#' are in 'symbol' format. It can be SummarizedExperiment, ExpressionSet, or matrix objects.
#' @param PCAmodel PCAGenomicSignatures object. Output from \code{buildAvgLoading}
#' function, a matrix of avaerage loadings, can be directly provided.
#' @return A list containing the score matrices for input datasets. Scores are
#' assigned to each sample (row) for each cluster (column).
#'
#' @export
calculateScore <- function(dataset, PCAmodel) {

    if (is(PCAmodel, "PCAGenomicSignatures")) {
        avg.loadings <- assay(PCAmodel)
    } else if (class(PCAmodel) %in% c("data.frame", "matrix")) {
        avg.loadings <- PCAmodel
    }

    if (!is.list(dataset)) {dataset <- list(dataset)}
    lapply(dataset, function(dat) {
        if (is(dat, "ExpressionSet")) {
            count <- Biobase::exprs(dat)
        } else if (is(dat, "SummarizedExperiment")) {
            count <- SummarizedExperiment::assay(dat)
        } else if (is(dat, "matrix")) {
            count <- dat
        }

        count <- count[apply(count, 1, function(x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]
        count <- apply(count, 1, function(x) {x - mean(x)}) %>% t
        gene_common <- intersect(rownames(avg.loadings), rownames(count))

        score <- t(count[gene_common,]) %*% apply(avg.loadings[gene_common,], 2,
                                                 function(x) x / sqrt(sum(x^2, na.rm = TRUE)))
        ## CRC paper version
        # score <- t(count[gene_common,]) %*% as.matrix(avg.loadings[gene_common,])
        # score <- (t(score) / apply(score, 2, sd)) %>% t

        colnames(score) <- colnames(avg.loadings)
        return(score)
    })
}
