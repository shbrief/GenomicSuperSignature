#' Find the studies contributing each RAV
#'
#' @param RAVmodel PCAGenomicSignatures object.
#' @param ind A numeric vector containing the RAV indexes. Under the default
#' (\code{NULL}), studies associated with all the RAV indexes will be returned
#' as a list.
#'
#' @return A list of character vectors. Under the default condition
#' (\code{ind = NULL}), all the RAVs will be checked for their contributing
#' studies and the length of the list will be same as the number of RAVs
#' (= \code{metadata(x)$k}). If you provide the \code{ind} argument, studies
#' associated with only the specified RAVs will be returned.
#'
#' @examples
#' data(miniRAVmodel)
#' findStudiesInCluster(miniRAVmodel, 1076)
#'
#' @note Mainly used for model building, within \link{buildAvgLoading}.
#' @export
findStudiesInCluster <- function(RAVmodel, ind = NULL) {

  if (!is.null(ind)) {.availableRAV(RAVmodel, ind)}
  
  if (is(RAVmodel,"PCAGenomicSignatures")) {
    x <- S4Vectors::metadata(RAVmodel)
    k <- x$k   # the number of clusters
  } else if (is.list(RAVmodel)) {   # [internal] this is for model building
    x <- RAVmodel
    x$size <- table(RAVmodel$cluster)
    k <- length(unique(RAVmodel$cluster))
  }

  # z is a binary matrix showing the cluster membership of PCs
  z <- matrix(0, ncol = k, nrow = sum(x$size))
  for (i in seq_along(x$cluster)) {
    z[i, x$cluster[i]] <- 1
  }
  colnames(z) <- paste0("Cl", k, "_",
                        formatC(seq_len(k), width=2, format="d", flag="0"))
  rownames(z) <- names(x$cluster)

  if (is.null(ind)) {
    studies <- list()
    for (i in seq_len(ncol(z))) {
      studies[[i]] <- rownames(z)[which(z[,i] == 1)]
      studies[[i]] <- lapply(studies[[i]],
                             function(x) {unlist(strsplit(x, "\\."))[1]})
      studies[[i]] <- unique(unlist(studies[[i]]))
      names(studies)[i] <- colnames(z)[i]
      res <- studies
    }
  } else {
    for (i in ind) {
      studies <- rownames(z[which(z[,i] == 1),])
      studies <- lapply(studies,
                        function(x) {unlist(strsplit(x, "\\."))[1]})
      res <- studies <- unique(unlist(studies))
    }
  }
  return(res)
}
