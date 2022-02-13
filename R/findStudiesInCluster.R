#' Find the studies contributing each RAV
#'
#' @param RAVmodel PCAGenomicSignatures object.
#' @param ind A numeric vector containing the RAV indexes. Under the default
#' (\code{NULL}), studies associated with all the RAV indexes will be returned
#' as a list.
#' @param studyTitle Default is \code{FALSE}. This parameter is effective only when
#' the \code{index} value is specificed. If it's \code{TRUE}, the output will be
#' a data frame with the study
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
findStudiesInCluster <- function(RAVmodel, ind = NULL, studyTitle = FALSE) {

  ## Check the input validity
  if (!is.null(ind)) {.availableRAV(RAVmodel, ind)}

  ## Extract cluster information: size and the number of PCs
  if (is(RAVmodel,"PCAGenomicSignatures")) {
    x <- S4Vectors::metadata(RAVmodel)
    k <- x$k   # the number of clusters
  } else if (is.list(RAVmodel)) {   # [internal] this is for model building
    x <- RAVmodel
    x$size <- table(RAVmodel$cluster)
    k <- length(unique(RAVmodel$cluster))
  }

  ## z is the binary matrix showing the cluster membership of PCs
  z <- matrix(0, ncol = k, nrow = sum(x$size))
  for (i in seq_along(x$cluster)) {
    z[i, x$cluster[i]] <- 1
  }
  colnames(z) <- paste0("Cl", k, "_",
                        formatC(seq_len(k), width=2, format="d", flag="0"))
  rownames(z) <- names(x$cluster)

  ## Extract study accession number part from PC names in each cluster
  if (is.null(ind)) {
    studies <- list()
    for (i in seq_len(ncol(z))) {
      studies[[i]] <- rownames(z)[which(z[,i] == 1)]
      studies[[i]] <- lapply(studies[[i]],
                             function(x) {unlist(strsplit(x, "\\.PC"))[1]})
      studies[[i]] <- unique(unlist(studies[[i]]))  # remove redundancy
      names(studies)[i] <- colnames(z)[i]
    }
  } else {
    studies <- .varByPCsInCluster(RAVmodel, ind)
  }

  ## Extract study titles
  if (!studyTitle) {
    res <- studies
  } else {
    studyMeta <- .getStudyMeta(RAVmodel)
    studyMeta <- studyMeta[,c("studyName", "title")]
    res <- dplyr::left_join(studies, studyMeta, by = "studyName")
  }

  return(res)
}
