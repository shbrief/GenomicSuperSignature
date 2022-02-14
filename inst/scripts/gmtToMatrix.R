#' Convert binary membership matrix into list
#'
#' @param x A binary membership matrix. Values should be 0 or 1.
#' @return A named list with the length of input matrix's columns, where each
#' element's names is same as column name of the input matrix.
#'
#' @export
matrixToList <- function(x) {
    res <- vector(mode = "list", length = ncol(x))
    names(res) <- colnames(x)
    for (i in seq_len(ncol(x))) {
        ind <- which(x[,i] == 1)
        res[[i]] <- rownames(x)[ind]
    }
    return(res)
}

#' Convert binary membership matrix into list
#'
#' @param x A binary membership matrix. Values should be 0 or 1.
#' @return A data frame with 2 columns: gs_name and entrez_gene
#'
#' @export
matrixToTERM2GENE <- function(x) {
    res <- as.data.frame(matrix(ncol = 2, nrow = 0))
    colnames(res) <- c("gs_name", "entrez_gene")

    for (i in seq_len(ncol(x))) {
        ind <- which(x[,i] == 1)
        res_sub <- data.frame(gs_name = colnames(x)[i],
                              entrez_gene = rownames(x)[ind])
        res <- rbind(res, res_sub)
    }
    return(res)
}




#' Combine all the prior knowledge matrices into one matrix
#'
#' @param priors A list of prior knowledge matrices.
#' @return A matrix of all the provided genesets. Row represents genes and column
#' represents pathways.
#'
#' @export
combinePriors <- function(priors) {

    genes <- character()
    for (i in seq_along(priors)) {
        genes <- c(genes, rownames(priors[[i]]))
        priors[[i]] <- as.data.frame(priors[[i]])
    }

    genes <- unique(genes)
    mat <- matrix(nrow = length(genes), ncol = 0)
    rownames(mat) <- genes

    for (i in seq_along(priors)) {
        priors[[i]] <- as.matrix(priors[[i]][genes, ])
        mat <- cbind(mat, priors[[i]])
    }

    mat[is.na(mat)] <- 0
    return(mat)
}


#' Generate prior knowledge matrix from gmt
#'
#' This function converts \code{.gmt} files into binary matrix.
#'
#' @importFrom EnrichmentBrowser getGenesets
#' @importFrom safe getCmatrix
#'
#' @param gmt.file A list of GMT files to build a prior knowledge matrix.
#' @param min.genes A minimum number of genes per genesets to be considered. Default is 10.
#' @param max.genes A maximum number of genes per genesets to be considered. Default is 500.
#'
#' @return A binary matrix with the geneset membership. Row represents genes and
#' columne represents pathways.
#'
#' @export
gmtToMatrix <- function(gmt.file, min.genes = 10, max.genes = 500) {

  genesets <- list()
  for (gmt.file in gmt.file) {
    gmtfile <- getGenesets(gmt.file)
    gmtMatrix <- getCmatrix(gmtfile, as.matrix = TRUE)
    genesetName <- gsub(".gmt", "", basename(gmt.file))
    genesets[[genesetName]] <- gmtMatrix
  }

  allPaths <- combinePriors(genesets)   # need to rewrite PLIER function as GSS function
  ind.exclude <- which(colSums(allPaths != 0) > max.genes | colSums(allPaths != 0) < min.genes)
  if (length(ind.exclude) != 0) {
    allPaths <- allPaths[, -ind.exclude]
  }
  return(allPaths)

}
