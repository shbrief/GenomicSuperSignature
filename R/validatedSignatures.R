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



#' Validation result in data frame
#'
#' @param data An output matrix from \code{\link{validate}} function.
#' @param num.out A number of highly validated PCclusters to output. Default is 5.
#' If any of the cutoff parameters are provided, \code{num.out} or the number of
#' filtered PCclusters, whichever smaller, will be chosen.
#' @param scoreCutoff A numeric value for the minimum correlation.
#' @param swCutoff A numeric value for the minimum average silhouette width.
#' @param clsizeCuoff A integer value for the minimum cluster size.
#'
#' @return A subset of the input matrix, which meets the given condition.
#' @export
validatedSignatures <- function(data, num.out = 5, scoreCutoff = NULL, swCutoff = NULL, clsizeCutoff = NULL) {
  ind <- which(rownames(data) == "score")

  if (!is.null(scoreCutoff)) {score_subset <- .filterByScore(data, scoreCutoff)} else {score_subset <- data}
  if (!is.null(swCutoff)) {sw_subset <- .filterByAvgSW(data, swCutoff)} else {sw_subset <- data}
  if (!is.null(clsizeCutoff)) {clsize_subset <- .filterByClusterSize(data, clsizeCutoff)} else {clsize_subset <- data}

  common_col <- intersect(colnames(score_subset),colnames(sw_subset)) %>%
    intersect(., colnames(clsize_subset))
  common_subset <- data[, common_col]
  ordered_ind <- order(common_subset[ind,], decreasing = TRUE)[1:min(num.out, ncol(common_subset))]
  dat <- common_subset[,ordered_ind]
  return(dat)
}
