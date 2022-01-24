#' Subset enriched pathways of RAV
#'
#' @import methods
#'
#' @param RAVmodel PCAGenomicSignatures object. Also an output from
#' \code{\link[clusterProfiler]{GSEA}} can be used.
#' @param ind A numeric vector containing the RAV number you want to check
#' enriched pathways. If not specified, this function returns results from all
#' the RAVs.
#' @param n The number of top and bottom pathways to be selected based on
#' normalized enrichment score (NES).
#' @param both Default is \code{FALSE}, where only the top \code{n} pathways
#' will be printed. If it is set to \code{TRUE}, the output will contain both
#' top and bottom \code{n} pathways.
#' @param include_nes Defalt is \code{FALSE}. If it set to \code{TRUE}, the
#' output will include both description and NES of the enriched pathway.
#'
#' @return A DataFrame with top and bottom \code{n} pathways from the
#' enrichment results.
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
subsetEnrichedPathways <- function(RAVmodel, ind = NULL, n = 10,
                                   both = FALSE, include_nes = FALSE) {

  ## Take the GSEA annotation
  if (is(RAVmodel, "PCAGenomicSignatures")) {
    gsea_loading <- gsea(RAVmodel)
  } else if(is.data.frame(RAVmodel) | is.matrix(RAVmodel) | is.list(RAVmodel)) {
    gsea_loading <- RAVmodel
  }

  ## RAVs to subset GSEA annotation from
  res <- list()
  if (is.null(ind)) {
    names <- names(gsea_loading)
  } else {names <- paste0("RAV", ind)}

  ## Subset the enriched pathways with top/bottom NES
  for (name in names) {
    x <- gsea_loading[[name]]
    if (is.null(dim(x))) {
      res[[name]] <- "There is no enriched pathway"
    } else {
      y <- x[,c("Description", "NES")]
      up <- y[order(y$NES, decreasing=TRUE),][seq_len(n),]
      rownames(up) <- paste0("Up_", seq_len(n))
      down <- y[order(y$NES, decreasing=FALSE),][seq_len(n),]
      rownames(down) <- paste0("Down_", seq_len(n))
      res[[name]] <- rbind(up, down)
    }
  }

  ## Remove NES if it is requested
  if (!include_nes) {
    for (name in names) {
      res[[name]] <- res[[name]][,"Description", drop = FALSE]
      colnames(res[[name]]) <- paste(name, colnames(res[[name]]), sep = ".")
    }
  }

  ## Combine list into one data frame
  res <- data.frame(res, check.names = FALSE)
  if (!both) {res <- res[seq_len(n),,drop=FALSE]}

  res <- S4Vectors::DataFrame(res)
  return(res)
}
