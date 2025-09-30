.set_n <- function(dataset, n) {
    if (is.null(n)) {
      if (ncol(dataset) >= 8) {
        n = 8
      } else {
        n = ceiling(ncol(dataset)*0.5)
      }
    }
    return(n)
}


#' Validating new dataset
#'
#' @importFrom irlba irlba
#' 
#' @param dataset A gene expression profile to be validated. Different classes
#' of objects can be used including ExpressionSet, SummarizedExperiment,
#' RangedSummarizedExperiment, or matrix. Rownames (genes) should be in 
#' human gene symbol format (HGNC). If dataset is a matrix, genes should be in 
#' rows and samples in columns. RNA-seq counts should be log(count + 1) 
#' prior to the `validate()` call.
#' @param avgLoading A matrix with genes by RAVs.
#' @param n A integer. The number of PCs to use for validation. It should be
#' equal or less then the number of samples in the input dataset. 
#' @param method A character string indicating which correlation coefficient is
#' to be computed. One of "pearson" (default), "kendall", or "spearman": can be
#' abbreviated.
#' @param scale Default is \code{FALSE}. If it is set to \code{TRUE}, rows 
#' will be converted to z-score prior to PCA. 
#'
#' @return A matrix of Pearson correlation coefficient (default, defined through
#' \code{method} argument) between RAVs (row) and the top 8 PCs from the
#' datasets (column)
#' 
#' @keywords internal
#'
.loadingCor <- function(dataset, 
                        avgLoading,
                        n = n,
                        method = method, 
                        scale = FALSE) {

    # Extract expression matrix from different classes
    dat <- .extractExprsMatrix(dataset)

    # # row normalization # Plan to remove `scale` argument <<<<<<<<<<<<<<<<<<<< 
    # stopifnot(length(scale) == 1L, !is.na(scale), is.logical(scale))
    # if (scale) {dat <- t(scale(t(dat)))}

    if (!is(dataset, "SingleCellExperiment")) {
        dat <- dat[apply(dat, 1,
                         function (x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]
        gene_common <- intersect(rownames(avgLoading), rownames(dat))
        prcomRes <- stats::prcomp(t(dat[gene_common,]))  # centered, but not scaled
        loadings <- prcomRes$rotation[, seq_len(n)]
    } else {    # Dimensional reduction for scRNAseq data
        gene_common <- intersect(rownames(avgLoading), rownames(dat))
        res <- irlba::irlba(as.matrix(t(dat[gene_common,])), 
                            nv = n, 
                            center = TRUE, scale = FALSE)
        loadings <- res$v
        rownames(loadings) <- rownames(dat[gene_common,])
        colnames(loadings) <- paste0("PC", seq_len(n))
    }
    
    loading_cor <- abs(stats::cor(avgLoading[gene_common,],
                                  loadings[gene_common,],
                                  use = "pairwise.complete.obs",
                                  method = method))
    return(loading_cor)
}


#' Validate new datasets
#'
#' @param dataset Single or a named list of SummarizedExperiment
#' (RangedSummarizedExperiment, ExpressionSet or matrix) object(s). Columns 
#' should contain samples (or cells for single-cell data), and the row names
#' should be in 'gene symbol' format. 
#' @param RAVmodel PCAGenomicSignatures object.
#' @param method A character string indicating which correlation coefficient is
#' to be computed. One of "pearson" (default), "kendall", or "spearman": can be
#' abbreviated.
#' @param maxFrom Select whether to display the maximum value from dataset's PCs
#' or avgLoadings. Under the default (\code{maxFrom="PC"}), the maximum
#' correlation coefficient from top 8 PCs for each avgLoading will be selected
#' as an output. If you choose (\code{maxFrom="avgLoading"}), the avgLoading
#' with the maximum correlation coefficient with each PC will be in the output.
#' @param n A integer. The number of PCs to use for validation. It should be
#' equal or less then the number of samples in the input dataset. If there are 
#' >= 8 samples in the input dataset or the input datasets is a list, it is set 
#' to 8 as a default. If there are less than 8 samples in the input, it is set 
#' to half of the number of samples as a default.
#' @param level Output format of validated result. Two options are available:
#' \code{c("max", "all")}. Default is "max", which outputs the matrix containing
#' only the maximum coefficient. To get the coefficient of all PCs, set this
#' argument as "all". \code{level = "all"} can be used only for one dataset.
#' @param scale Default is \code{FALSE}. If it is set to \code{TRUE}, dataset
#' will be row normalized.
#'
#' @return A data frame containing the maximum pearson correlation coefficient
#' between the top `n` PCs of the dataset and pre-calculated average loadings
#' (in row) of training datasets (\code{score} column). It also contains other
#' metadata associated with each RAV: \code{PC} for one of the top `n` PCs of 
#' the dataset that results in the given \code{score}, \code{sw} for the average
#' silhouette width of the RAV, \code{cl_size} for the size of each RAV.
#'
#' If the input for \code{dataset} argument is a list of different datasets,
#' each row of the output represents a new dataset for test, and each column
#' represents clusters from training datasets. If \code{level = "all"}, a list
#' containing the matrices of the pearson correlation coefficient between all
#' top 8 PCs of the datasets and avgLoading.
#'
#' @examples
#' data(miniRAVmodel)
#' library(bcellViper)
#' data(bcellViper)
#' validate(dset, miniRAVmodel)
#' validate(dset, miniRAVmodel, maxFrom = "avgLoading")
#'
#' @export
validate <- function(dataset, 
                     RAVmodel, 
                     n = NULL,
                     method = "pearson",
                     maxFrom = "PC", 
                     level = "max", 
                     scale = FALSE) {
  
    if (!is.list(dataset)) {
        n <- .set_n(dataset, n)
        if (ncol(dataset) < n) {
            stop("n should be equal or less than the number of samples.")}
    } else {
        if (is.null(n)) {n <- 8}
        if (any(lapply(dataset, ncol) < n)) {
            stop("n should be equal or less than the number of samples.")}
        if (level == "all") {
            stop("'level = \"all\"' is not available for a list of datasets.")}
    }

    sw <- silhouetteWidth(RAVmodel)
    cl_size <- S4Vectors::metadata(RAVmodel)$size
    avgLoading <- SummarizedExperiment::assay(RAVmodel)

    # The maximum correlation coefficient among PCs
    if (maxFrom == "PC") {
        # For a single dataset
        if (!is.list(dataset)) {
            x <- .loadingCor(dataset, avgLoading, n, method, scale)
            if (level == "max") {
                z <- apply(x, 1, max) %>% as.data.frame   # rowMax
                z$PC <- apply(x, 1, which.max)
                colnames(z)[1] <- "score"
                z$sw <- sw   # Silhouette width
                z$cl_size <- cl_size   # Cluster size
                z$cl_num <- readr::parse_number(rownames(z))   # Cluster number
                res <- z
                return(res)
            } else if (level == "all") {
                res <- x
                return(res)
            }
        } else {
            # For a list of datasets
            x <- lapply(dataset, .loadingCor, avgLoading, n, method, scale)
            l <- nrow(x[[1]]) # the number of RAVs in validation output
            if (level == "max") {
                z <- vapply(x, function(y) {apply(y, 1, max)},
                            FUN.VALUE = numeric(l))
                zPC <- vapply(x, function(y) {apply(y, 1, which.max)},
                              FUN.VALUE = integer(l))
                colnames(zPC) <- paste0(colnames(zPC), "_PC")
                res <- cbind(z, zPC)
                return(res)
            } else if (level == "all") {
                res <- x
                return(res)
            }
        }
    }

    # The maximum correlation coefficient among avgLoadings
    if (maxFrom == "avgLoading") {
        if (!is.list(dataset)) {
            x <- .loadingCor(dataset, avgLoading, n, method)
            if (level == "max") {
                z <- apply(x, 2, max) %>% as.data.frame # colMax
                max_z_ind <- apply(x, 2, which.max)
                z$validated_RAV <- rownames(x)[max_z_ind]
                colnames(z)[1] <- "score"
                return(z)
            } else if (level == "all") {
                return(x)
            }
        } else {
            x <- lapply(dataset, .loadingCor, avgLoading, method)
            if (level == "max") {
                z <- apply(x, 2, max)
                return(z)
            } else if (level == "all") {
                return(x)
            }
        }
    }
}
