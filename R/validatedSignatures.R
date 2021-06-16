.filterBy <- function(data, cutoff, type = c("cl_size", "score", "sw")) {
  type <- match.arg(type)
  ind <- rownames(data) == type
  val_ind <- which(data[ind,] >= cutoff)
  data[,val_ind]
}

.validatedSignaturesForMulipleStudies <- function(data, num.out = num.out,
                                                  scoreCutoff = NULL) {
  studies <- rownames(data)
  ind <- grep("_PC", studies)
  data <- data[-ind,]

  if (!is.null(scoreCutoff)) {scoreCutoff <- scoreCutoff}
  else {scoreCutoff <- 0.7}   # default cutoff for multiple-studies case

  above_cutoff <- apply(data, 2, function(x) {any(x > scoreCutoff)})
  data <- data[, above_cutoff, drop=FALSE]
  return(data)
}


#' Validation result in data frame
#'
#' @param val_all An output matrix from \code{\link{validate}} function. If this
#' input is from multiple datasets, only \code{scoreCutoff} argument will be
#' considered and other inputs will be ignored.
#' @param num.out A number of highly validated RAVs to output. Default is 5.
#' If any of the cutoff parameters are provided, \code{num.out} or the number of
#' filtered RAVs, whichever smaller, will be chosen.
#' @param scoreCutoff A numeric value for the minimum correlation. For
#' multi-studies case, the default is 0.7.
#' @param swCutoff A numeric value for the minimum average silhouette width.
#' @param clsizeCutoff An integer value for the minimum cluster size.
#' @param indexOnly A logical. Under the default (= \code{FALSE}), the detailed
#' information on validated RAVs, such as score, average silhouette width,
#' cluster size, is printed. If it is set TRUE, only the RAV number will be
#' printed.
#' @param whichPC An integer value between 1 and 8. PC number of your data to
#' check the validated signatures with. Under the default (\code{NULL}), it
#' outputs top scored signatures with any PC of your data.
#'
#' @return A subset of the input matrix, which meets the given condition.
#'
#' @examples
#' data(miniRAVmodel)
#' library(bcellViper)
#' data(bcellViper)
#' val_all <- validate(dset, miniRAVmodel)
#' validatedSignatures(val_all, num.out = 3, scoreCutoff = 0)
#'
#' @export
validatedSignatures <- function(val_all, num.out = 5, scoreCutoff = NULL,
                                swCutoff = NULL, clsizeCutoff = NULL,
                                indexOnly = FALSE, whichPC = NULL) {
  # Input validation
  stopifnot(length(indexOnly) == 1L, !is.na(indexOnly), is.logical(indexOnly))

  data <- t(val_all)
  ind <- which(rownames(data) == "score")
  pc_ind <- which(rownames(data) == "PC")

  # If the validation result is from the list of datasets
  if (length(ind) == 0) {
    res <- .validatedSignaturesForMulipleStudies(data, scoreCutoff=scoreCutoff)

    if (!indexOnly) {
      return(res)
    } else {
      res_ind <- gsub("RAV", "", colnames(res))
      res_ind <- as.numeric(res_ind)
      return(res_ind)
    }
    stop
  }

  if (!is.null(whichPC)) {
    if (!whichPC %in% seq_len(8)) {
      stop("whichPC should be an integer between 1 and 8.")
    }
    subset_ind <- which(data[pc_ind,] == whichPC)
    data <- data[,subset_ind, drop = FALSE]
  }

  if (!is.null(scoreCutoff)) {
    score_subset <- .filterBy(data, scoreCutoff, "score")
  } else {score_subset <- data}
  if (!is.null(swCutoff)) {
    sw_subset <- .filterBy(data, swCutoff, "sw")
  } else {sw_subset <- data}
  if (!is.null(clsizeCutoff)) {
    clsize_subset <- .filterBy(data, clsizeCutoff, "cl_size")
  } else {clsize_subset <- data}

  common_col <- intersect(colnames(score_subset), colnames(sw_subset))
  common_col <- intersect(common_col, colnames(clsize_subset))

  if (length(common_col) == 0) {
    dat <- data[, 0]
  } else {
    common_subset <- data[, common_col]
    ordered_ind <- order(common_subset[ind,], decreasing = TRUE)
    n_out <- min(num.out, ncol(common_subset))
    ordered_ind <- ordered_ind[seq_len(n_out)]
    dat <- common_subset[,ordered_ind, drop = FALSE]
  }

  if (!indexOnly) {
    return(t(dat))
  } else {
    validatedIndex <- dat["cl_num",]
    validatedIndex <- as.numeric(validatedIndex)
    return(validatedIndex)
  }
}
