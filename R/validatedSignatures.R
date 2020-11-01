.filterByClusterSize <- function(data, cutoff) {
  ind <- which(rownames(data) == "cl_size")
  val_ind <- which(data[ind,] >= cutoff)
  dat <- data[,val_ind]
  return(dat)
}

.filterByScore <- function(data, cutoff) {
  ind <- which(rownames(data) == "score")
  val_ind <- which(data[ind,] >= cutoff)
  dat <- data[,val_ind]
  return(dat)
}

.filterByAvgSW <- function(data, cutoff) {
  ind <- which(rownames(data) == "sw")
  val_ind <- which(data[ind,] >= cutoff)
  dat <- data[,val_ind]
  return(dat)
}

.validatedSignaturesForMulipleStudies <- function(data, num.out = num.out, scoreCutoff = NULL) {
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
#' input is from multiple datasets, only \code{scoreCutoff} argument will be considered
#' and other inputs will be ignored.
#' @param num.out A number of highly validated RAVs to output. Default is 5.
#' If any of the cutoff parameters are provided, \code{num.out} or the number of
#' filtered RAVs, whichever smaller, will be chosen.
#' @param scoreCutoff A numeric value for the minimum correlation. For multi-studies
#' case, the default is 0.7.
#' @param swCutoff A numeric value for the minimum average silhouette width.
#' @param clsizeCutoff An integer value for the minimum cluster size.
#' @param indexOnly A logical. Under the default (= FALSE), the detailed information
#' on validated RAVs, such as score, average silhouette width, cluster size, is
#' printed. If it is set TRUE, only the RAV number will be printed.
#' @param whichPC An integer value between 1 and 8. PC number of your data to check
#' the validated signatures with. Under the default (\code{NULL}), it outputs top
#' scored signatures with any PC of your data.
#'
#' @return A subset of the input matrix, which meets the given condition.
#'
#' @examples
#' data(miniPCAmodel)
#' library(bcellViper)
#' data(bcellViper)
#' val_all <- validate(dset, miniPCAmodel)
#' validatedSignatures(val_all, num.out = 3, scoreCutoff = 0)
#' #             score PC          sw cl_size cl_num
#' # RAV1076 0.5950767  2 -0.04447124      10      1
#' # RAV2538 0.5838616  2  0.06996166       4      2
#' # RAV338  0.5709072  2 -0.04683319      21      3
#'
#' @export
validatedSignatures <- function(val_all, num.out = 5, scoreCutoff = NULL, swCutoff = NULL,
                                clsizeCutoff = NULL, indexOnly = FALSE, whichPC = NULL) {
  data <- t(val_all)
  ind <- which(rownames(data) == "score")
  pc_ind <- which(rownames(data) == "PC")

  # If the validation result is from the list of datasets
  if (length(ind) == 0) {
    res <- .validatedSignaturesForMulipleStudies(data, scoreCutoff = scoreCutoff)

    if (isFALSE(indexOnly)) {
      return(res)
    } else {
      res_ind <- gsub("RAV", "", colnames(res))
      res_ind <- as.numeric(res_ind)
      return(res_ind)
    }
    stop
  }

  if (!is.null(whichPC)) {
    if (!whichPC %in% c(1:8)) {stop("whichPC should be an integer between 1 and 8.")}
    subset_ind <- which(data[pc_ind,] == whichPC)
    data <- data[,subset_ind, drop = FALSE]
  }

  if (!is.null(scoreCutoff)) {score_subset <- .filterByScore(data, scoreCutoff)} else {score_subset <- data}
  if (!is.null(swCutoff)) {sw_subset <- .filterByAvgSW(data, swCutoff)} else {sw_subset <- data}
  if (!is.null(clsizeCutoff)) {clsize_subset <- .filterByClusterSize(data, clsizeCutoff)} else {clsize_subset <- data}

  common_col <- intersect(colnames(score_subset), colnames(sw_subset))
  common_col <- intersect(common_col, colnames(clsize_subset))

  if (length(common_col) == 0) {
    dat <- data[, 0]
  } else {
    common_subset <- data[, common_col]
    ordered_ind <- order(common_subset[ind,], decreasing = TRUE)[1:min(num.out, ncol(common_subset))]
    dat <- common_subset[,ordered_ind, drop = FALSE]
  }

  if (isFALSE(indexOnly)) {
    return(t(dat))
  } else {
    validatedIndex <- dat["cl_num",]
    validatedIndex <- as.numeric(validatedIndex)
    return(validatedIndex)
  }
}
