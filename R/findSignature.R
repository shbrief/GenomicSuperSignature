#' Find the PCclusters with the keyword-containing enriched pathways
#'
#' This function finds PCclusters containing the keyword you provide. If you provide
#' "the number of keyword-containing pathways per PCcluster" in argument \code{k},
#' it will give you the PCcluster number.
#'
#' @param PCAmodel PCAGenomicSignatures-object
#' @param keyword A character vector. If you are searching for multiple keywords
#' at the same time, use \code{\link{paste}} with \code{collapse = "|"} argument.
#' @param n The number of top ranked (based on abs(NES)) pathways you want to search your keyword
#' @param k The number of keyword-containing pathways you want to get the PCcluster
#' number. Under default (\code{NULL}), the output will be a data frame with two
#' columns: '# of keyword-containing pathways' and 'Freq'. If you assign the value
#' for this argument, the output will be an integer vector containing the PCcluster
#' index.
#'
#' @return A data frame or integer vector depending on the parameter \code{k}.
#'
#' @examples
#' # findSignature(PCAmodel, "CD4|T cell")
#' # findSignature(PCAmodel, "CD4|T cell", k = 2)
#'
#' @export
findSignature <- function(PCAmodel, keyword, n = 5, k = NULL) {
  gsea_all <- gsea(PCAmodel)

  # top n pathways based on abs(NES) value
  topPathways <- sapply(gsea_all, function(x) {
    y <- x[order(abs(x$NES), decreasing = TRUE),][seq_len(n),]
    grep(keyword, y[["Description"]], ignore.case = TRUE)
  })

  nTopPathways <- sapply(topPathways, length)
  res <- nTopPathways %>% table %>% as.data.frame
  colnames(res)[1] <- c("# of keyword-containing pathways")

  if (is.null(k)) {
    return(res)
  } else {
    which(nTopPathways == k)
  }
}



#' Find the rank of your keyword in the PCcluster's GSEA annotation
#'
#' Once you provide PCAmodel, keyword you're searching for, and the PCcluster number
#' to this function, it will give you the abs(NES)-based rank of your keyword in
#' the enriched pathways of the target PCcluster. If can be useful to find out how
#' uniquely your keyword-containing pathways are represented.
#'
#' @param PCAmodel PCAGenomicSignatures-object.
#' @param keyword A character vector. If you are searching for multiple keywords
#' at the same time, use \code{\link{paste}} with \code{collapse = "|"} argument.
#' @param ind An integer. The PCcluster number you want to check.
#' @param n An interger. The number of top enriched pathways (based on abs(NES))
#' to search. Under default (\code{NULL}), all the enriched pathways are used.
#' @param includeTotal Default is \code{FALSE}. If it is set to \code{TRUE}, the
#' total number of enriched pathways will be also printed out as an output.
#'
#' @return A character containing the rank of keyword-containing pathways (separated
#' by |), followed by the total number of enriched pathways in parenthesis.
#'
findKeywordInPCcluster <- function(PCAmodel, keyword, ind, n = NULL, includeTotal = FALSE) {
  gsea <- gsea(PCAmodel)[[ind]]
  total <- nrow(gsea)

  if (is.null(n)) {
    gsea <- gsea[order(abs(gsea$NES), decreasing = TRUE),,drop = FALSE]
  } else {
    gsea <- gsea[order(abs(gsea$NES), decreasing = TRUE),,drop = FALSE][seq_len(n),]
  }

  keywordRank <- grep(keyword, gsea$Description, ignore.case = TRUE)
  if (isFALSE(includeTotal)) {
    res <- paste(keywordRank, collapse = "|")
  } else {
    res <- paste0(paste(keywordRank, collapse = "|"), " (out of ", total, ")")
  }
  return(res)
}





