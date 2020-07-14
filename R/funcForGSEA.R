#' @title Order genes in loading vectors
#' @description This function takes Z matrix (= average loadings) and orders the
#' genes in each loading vector (= PCcluster) in a descending manner.
#'
#' @param LoadingMatrix An avgloading matrix. Rows represent genes and columns
#' represent clusters of principle components
#' @param LoadingVector A list of column names or indexes of \code{LoadingMatrix}
#' you want to check. Default is \code{NULL}, under which the function takes the
#' all column names of \code{LoadingMatrix}
#' @param abs Under the defaul condition (\code{TRUE}), this function will
#' create a gene list based on the absolute value.
#'
#' @return A list of loadings selected by \code{LoadingVector}, where all the
#' genes in each loading are listed in descending order.
#'
#' @export
makeGeneList <- function(LoadingMatrix, LoadingVector = NULL, abs = TRUE) {

    if (is.null(LoadingVector)) {
        LoadingVector <- colnames(LoadingMatrix)
    }

    if (abs == TRUE) {
        LoadingMatrix <- abs(LoadingMatrix)
    }

    geneLists <- list()
    for (x in LoadingVector) {
        if (x %in% colnames(LoadingMatrix)) {
            geneList <- LoadingMatrix[, x]    # feature 1: numeric vector
            names(geneList) <- rownames(LoadingMatrix)    # feature 2: named vector
        } else {
            geneList <- LoadingMatrix[, c(x)]   # If index, not the name of column, is provided.
            names(geneList) <- rownames(LoadingMatrix)
        }

        geneList <- sort(geneList, decreasing = TRUE)    # feature 3: decreasing order
        geneLists[[as.character(x)]] <- geneList
    }
    return(geneLists)
}


#' @title GSEA on pre-ordered gene lists
#' @description This function is a wrapper of \code{\link[clusterProfiler]{GSEA}} function,
#' making it applicable to a list of gene lists. Set seed for reproducible result.
#'
#' @importFrom clusterProfiler GSEA
#'
#' @param geneList A list of genes ordered by rank
#' @param TERM2GENE User input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
#' @param TERM2NAME User input of TERM TO NAME mapping, a data.frame of 2 column with term and name. Optional.
#' @param minGSSize A mininum size of gene set to be analyzed
#' @param maxGSSize A maximum size of gene set to be analyzed
#' @param pvalueCutoff p-value cutoff
#' @param verbose Logical. Default is \code{FALSE}
#' @param ... Any additional argument inherited from \code{\link[clusterProfiler]{GSEA}}.
#'
#' @return A list of \code{gseaResult} objects. \code{NA} if there is no enrichment result.
#'
#' @export
run_gsea <- function(geneList, TERM2GENE, TERM2NAME,
                     minGSSize = 10, maxGSSize = 500,
                     pvalueCutoff = 0.05, verbose = FALSE, ...) {
    gsea <- list()
    for (x in names(geneList)) {
        res <- clusterProfiler::GSEA(geneList[[x]],
                                    TERM2GENE = TERM2GENE,
                                    TERM2NAME = TERM2NAME,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    verbose = verbose,
                                    pvalueCutoff = pvalueCutoff,
                                    ...)

        # collect GSEA outputs with enrichment result under the assigned pvalueCutoff
        if (nrow(res) != 0) {
            gsea[[x]] <- res
        }
    }
    return(gsea)
}


#' Subset enriched pathways of loading vectors
#'
#' This function is renamed from `topPathways` to `subsetPathways`.
#'
#' @importFrom enrichplot cnetplot
#' @import methods
#'
#' @param PCAmodel PCAGenomicSignatures object. Also an output from `clusterProfiler::GSEA` can be used.
#' @param ind A numeric vector containing the PCcluster number you want to check
#' enriched pathways. If not specified, this function returns results from all the PCclusters.
#' @param n The number of top and bottom pathways to be selected based on normalized
#' enrichment score (NES).
#' @param both Default is \code{FALSE}, where only the top \code{n} pathways will
#' be printed. If it is set to \code{TRUE}, the ouput will contain both top and
#' bottom \code{n} pathways.
#'
#'
#' @return A data frame with top and bottom \code{n} pathways from the enrichment results.
#'
#' @export
subsetPathways <- function(PCAmodel, ind = NULL, n = 10, both = FALSE) {

    if (is(PCAmodel, "PCAGenomicSignatures")) {
        gsea_loading <- unlist(gsea(PCAmodel))
    } else if (class(PCAmodel) %in% c("data.frame", "matrix", "list")) {
        gsea_loading <- PCAmodel
    }

    res <- list()
    for (name in names(gsea_loading)) {
        x <- gsea_loading[[name]]
        up <- x$Description[order(x$NES, decreasing=TRUE)][1:n]
        down <- x$Description[order(x$NES, decreasing=FALSE)][1:n]
        res[[name]] <- c(up, down)

        ### Adding NES and qvalues?
        # y <- x[,c("Description", "NES", "qvalues")]
        # up <- y[order(y$NES, decreasing=TRUE),][1:n,]
        # down <- y[order(y$NES, decreasing=FALSE),][1:n,]
        # res[[name]] <- rbind(up, down)
    }

    res <- data.frame(res, check.names = FALSE)
    rnames <- c(paste0("Up_", 1:n), paste0("Down_", 1:n))
    rownames(res) <- rnames

    if (both == FALSE) {res <- res[1:n,]}

    if (is.null(ind)) {
        res <- DataFrame(res)
        return(res)
    } else {
        col_num <- which(colnames(res) %in% c(paste0("PCcluster", ind)))
        res <- res[, col_num, drop = FALSE]
        res <- DataFrame(res)
        return(res)
    }
}
