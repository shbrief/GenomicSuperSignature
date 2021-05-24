#' Search the top enriched pathways for RAV
#'
#' @param RAVmodel PCAGenomicSignatures object.
#' @param ind An integer for RAV you want to check the enriched pathways.
#' @param n A number of top enriched pathways to output. Default is 5.
#' @param abs Default is \code{FALSE}. If it's set to \code{TRUE}, the enriched
#' pathways will be listed based on \code{abs(NES)}.
#'
#' @return A data frame with \code{n} rows and 4 columns;
#' Description, NES, pvalue, and qvalues
#'
#' @examples
#' data(miniRAVmodel)
#' annotateRAV(miniRAVmodel, ind = 695)
#'
#' @export
annotateRAV <- function(RAVmodel, ind, n = 5, abs = FALSE) {

  # extract GSEA results of the given cluster
  cl_name <- paste0("RAV", ind)
  annotatedCluster <- gsea(RAVmodel)[[cl_name]]

  # absolute value of NES
  FUN <- if (abs) {abs} else {I}

  # subset GSEA results
  topAnnotation <- annotatedCluster[
    order(FUN(annotatedCluster$NES), decreasing = TRUE),,drop = FALSE]
  topAnnotation <- topAnnotation[seq_len(n),]
  rownames(topAnnotation) <- NULL

  return(topAnnotation)
}
