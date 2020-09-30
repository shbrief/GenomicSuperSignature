#' Search the top enriched pathways for RAV
#'
#' @param PCAmodel PCAGenomicSignatures object.
#' @param ind An integer for RAV you want to check the enriched pathways.
#' @param n A number of top enriched pathways to output. Default is 5.
#' @param abs Default is \code{FALSE}. If it's set to \code{TRUE}, the enriched
#' pathways will be listed based on \code{abs(NES)}.
#'
#' @return A data frame with \code{n} rows and 4 columns; Description, NES, pvalue, and qval
#'
#' @export
annotatePCcluster <- function(PCAmodel, ind, n = 5, abs = FALSE) {
  cl_name <- paste0("PCcluster", ind)
  annotatedCluster <- gsea(PCAmodel)[[cl_name]]
  if (isTRUE(abs)) {
    topAnnotation <- annotatedCluster[order(abs(annotatedCluster$NES), decreasing = TRUE),,drop = FALSE][1:n,]
  } else {
    topAnnotation <- annotatedCluster[order(annotatedCluster$NES, decreasing = TRUE),,drop = FALSE][1:n,]
  }
  rownames(topAnnotation) <- NULL
  return(topAnnotation)
}
