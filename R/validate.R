#' Validating new dataset
#'
#' @param dataset A expression dataset to validate. Genes in rows and samples in
#' columns. Gene names should be in 'symbol' format. It can be ExpressionSet,
#' SummarizedExperiment, RangedSummarizedExperiment, or matrix.
#' @param avgLoading A matrix with genes by PCclusters.
#' @param method A character string indicating which correlation coefficient is
#' to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param scale Default is \code{FALSE}. If it is set to \code{TRUE}, dataset will
#' be row normalized by \link{rowNorm} function.
#'
#' @return A matrix of Pearson correlation coefficient (default, defined through \code{method}
#' argument) between PCclusters (row) and the top 8 PCs from the datasets (column)
#'
.loadingCor <- function(dataset, avgLoading, method = "pearson", scale = FALSE) {

    if (any(class(dataset) == "ExpressionSet")) {
        dat <- Biobase::exprs(dataset)
    } else if (any(class(dataset) %in% c("SummarizedExperiment", "RangedSummarizedExperiment"))) {
        dat <- SummarizedExperiment::assay(dataset)
    } else if (any(class(dataset) == "matrix")) {
        dat <- dataset
    } else {
        stop("'dataset' should be one of the following objects: ExpressionSet,
             SummarizedExperiment, RangedSummarizedExperiment, and matrix.")
    }

    if (isTRUE(scale)) {dat <- rowNorm(dat)}   # row normalization
    dat <- dat[apply(dat, 1, function (x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]
    gene_common <- intersect(rownames(avgLoading), rownames(dat))
    prcomRes <- stats::prcomp(t(dat[gene_common,]))  # centered, but not scaled by default
    loadings <- prcomRes$rotation[, 1:8]
    loading_cor <- abs(stats::cor(avgLoading[gene_common,], loadings[gene_common,],
                                  use = "pairwise.complete.obs",
                                  method = method))
    return(loading_cor)
}


#' Validate new datasets
#'
#' @importFrom dplyr tbl_df
#'
#' @param dataset Single or a list of SummarizedExperiment (RangedSummarizedExperiment,
#' ExpressionSet or matrix) object(s). Gene names should be in 'symbol' format.
#' @param PCAmodel PCAGenomicSignature object. You can also provide signature model matrix directly.
#' @param method A character string indicating which correlation coefficient is
#' to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param maxFrom Select whether to display the maximum value from dataset's PCs or avgLoadings.
#' Under the default (\code{maxFrom="PC"}), the maximum correlation coefficient from
#' top 8 PCs for each avgLoading will be selected as an output.If you choose (\code{maxFrom="avgLoading"}),
#' the avgLoading with the maximum correlation coefficient with each PC will be in the output.
#' @param level Defibe how to ouput validation in two different forms, \code{c("max", "all")}.
#' Default is "max", which outputs the matrix containing only the maximum coefficient.
#' To get the coefficient of all 8 PCs, set this argument as "all".
#' @param scale Default is \code{FALSE}. If it is set to \code{TRUE}, dataset will
#' be row normalized by \link{rowNorm} function.
#'
#' @return A data frame containing the maximum pearson correlation coefficient between
#' the top 8 PCs of the dataset and pre-calculated average loadings (in row) of training
#' datasets (\code{score} column). It also contains other metadata associated with
#' each PCcluster: \code{PC} for one of the top 8 PCs of the dataset that results
#' in the given \code{score}, \code{sw} for the average silhouette width of the PCcluster,
#' \code{cl_size} for the size of each PCcluster.
#' If the input for \code{dataset} argument is a list of different datasets, each row
#' of the output represents a new dataset for test, and each column represents
#' clusters from training datasets. If \code{level = "all"}, a list containing the matrices
#' of the pearson correlation coefficient between all top 8 PCs of the datasets and
#' avgLoading.
#'
#' @export
validate <- function(dataset, PCAmodel, method = "pearson",
                     maxFrom = "PC", level = "max", scale = FALSE)
    {
    sw <- silhouetteWidth(PCAmodel)
    cl_size <- S4Vectors::metadata(PCAmodel)$size

    if (class(PCAmodel) %in% c("PCAGenomicSignatures", "PLIERGenomicSignatures")) {
        avgLoading <- SummarizedExperiment::assay(PCAmodel)
    } else {avgLoading <- PCAmodel}

    # The maximum correlation coefficient among PCs
    if (maxFrom == "PC") {
        # For a single dataset
        if (!is.list(dataset)) {
            x <- .loadingCor(dataset, avgLoading, method, scale)
            if (level == "max") {
                z <- apply(x, 1, max) %>% as.data.frame   # rowMax
                z$PC <- apply(x, 1, which.max)
                colnames(z)[1] <- "score"
                z$sw <- sw   # Silhouette width
                z$cl_size <- cl_size   # Cluster size
                z$cl_num <- 1:nrow(z)   # Cluster number
                res <- z
            } else if (level == "all") {
                res <- x
            }
        } else {
            # For a list of datasets
            x <- lapply(dataset, .loadingCor, avgLoading, method, scale)
            if (level == "max") {
                z <- sapply(x, function(y) {apply(y, 1, max)})
                zPC <- sapply(x, function(y) {apply(y, 1, which.max)})
                colnames(zPC) <- paste0(colnames(zPC), "_PC")
                res <- cbind(z, zPC)
            } else if (level == "all") {
                res <- x
            }
        }
        # return(t(res))
        return(res)
    }

    # The maximum correlation coefficient among avgLoadings
    else if (maxFrom == "avgLoading") {
        if (!is.list(dataset)) {
            x <- .loadingCor(dataset, avgLoading, method)
            if (level == "max") {
                z <- apply(x, 2, max) %>% as.data.frame # colMax
                z$avgLoading <- apply(x, 2, which.max)
                colnames(z)[1] <- "score"
                return(t(z))
            } else if (level == "all") {
                return(t(x))
            }
        } else {
            x <- lapply(dataset, .loadingCor, avgLoading, method)
            if (level == "max") {
                z <- sapply(x, function(y) {apply(y, 2, max)})
                return(t(z))
            } else if (level == "all") {
                return(t(x))
            }
        }
    }
}
