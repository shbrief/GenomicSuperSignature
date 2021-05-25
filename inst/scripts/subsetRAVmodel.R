#' @param RAVmodel GenomicSignatures object
#' @param selectRAVs An integer vector. Index of RAVs to include.
#'
subsetRAVmodel <- function(RAVmodel, selectRAVs) {

  # Subset RAVindex
  miniRAVmodel <- RAVmodel[,selectRAVs]

  # Subset metadata$cluster
  cl_membership <- which(metadata(miniRAVmodel)$cluster %in% selectRAVs)
  metadata(miniRAVmodel)$cluster <-metadata(miniRAVmodel)$cluster[cl_membership]

  # Subset metadata$size
  metadata(miniRAVmodel)$size <- metadata(miniRAVmodel)$size[ind_all]
}
