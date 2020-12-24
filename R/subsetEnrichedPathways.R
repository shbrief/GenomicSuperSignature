#' Subset enriched pathways of RAV
#'
#' This function finds the top enriched pathways of a given RAV. This function is
#' renamed from \code{subsetPathways} to \code{subsetEnrichedPathways}.
#'
#' @import methods
#'
#' @param RAVmodel PCAGenomicSignatures object. Also an output from \code{\link[clusterProfiler]{GSEA}} can be used.
#' @param ind A numeric vector containing the RAV number you want to check
#' enriched pathways. If not specified, this function returns results from all the RAVs.
#' @param n The number of top and bottom pathways to be selected based on normalized
#' enrichment score (NES).
#' @param both Default is \code{FALSE}, where only the top \code{n} pathways will
#' be printed. If it is set to \code{TRUE}, the ouput will contain both top and
#' bottom \code{n} pathways.
#'
#' @return A DataFrame with top and bottom \code{n} pathways from the enrichment results.
#'
#' @export
subsetEnrichedPathways <- function(RAVmodel, ind = NULL, n = 10, both = FALSE) {

  if (is(RAVmodel, "PCAGenomicSignatures")) {
    gsea_loading <- gsea(RAVmodel)
  } else if (class(RAVmodel) %in% c("data.frame", "matrix", "list")) {
    gsea_loading <- RAVmodel
  }

  res <- list()
  for (name in names(gsea_loading)) {
    x <- gsea_loading[[name]]
    if (is.null(dim(x))) {
      res[[name]] <- "There is no enriched pathway"
    } else {
      up <- x$Description[order(x$NES, decreasing=TRUE)][1:n]
      down <- x$Description[order(x$NES, decreasing=FALSE)][1:n]
      res[[name]] <- c(up, down)
    }

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
    res <- S4Vectors::DataFrame(res)
    return(res)
  } else {
    col_num <- which(colnames(res) %in% c(paste0("RAV", ind)))
    res <- res[, col_num, drop = FALSE]
    res <- S4Vectors::DataFrame(res)
    return(res)
  }
}
