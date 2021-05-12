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
#' be printed. If it is set to \code{TRUE}, the output will contain both top and
#' bottom \code{n} pathways.
#'
#' @return A DataFrame with top and bottom \code{n} pathways from the enrichment results.
#'
#' @examples
#' data(miniRAVmodel)
#'
#' # all RAVS in model
#' subsetEnrichedPathways(miniRAVmodel,n=5)
#'
#' # only a specific RAV (note the colnames above)
#' subsetEnrichedPathways(miniRAVmodel,ind=695,n=5)
#'
#'
#' @export
subsetEnrichedPathways <- function(RAVmodel, ind = NULL, n = 10, both = FALSE) {

  if (is(RAVmodel, "PCAGenomicSignatures")) {
    gsea_loading <- gsea(RAVmodel)
  } else if (is.data.frame(RAVmodel) | is.matrix(RAVmodel) | is.list(RAVmodel)) {
    gsea_loading <- RAVmodel
  }

  res <- list()
  for (name in names(gsea_loading)) {
    x <- gsea_loading[[name]]
    if (is.null(dim(x))) {
      res[[name]] <- "There is no enriched pathway"
    } else {
      up <- x$Description[order(x$NES, decreasing=TRUE)][seq_len(n)]
      down <- x$Description[order(x$NES, decreasing=FALSE)][seq_len(n)]
      res[[name]] <- c(up, down)
    }

    ### Adding NES and qvalues?
    # y <- x[,c("Description", "NES", "qvalues")]
    # up <- y[order(y$NES, decreasing=TRUE),][seq_len(n),]
    # down <- y[order(y$NES, decreasing=FALSE),][seq_len(n),]
    # res[[name]] <- rbind(up, down)
  }

  res <- data.frame(res, check.names = FALSE)
  rnames <- c(paste0("Up_", seq_len(n)), paste0("Down_", seq_len(n)))
  rownames(res) <- rnames

  if (!both) {res <- res[seq_len(n),]}

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
