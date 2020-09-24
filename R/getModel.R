#' Download the PCAmodel from GCP
#'
#' @importFrom AnVIL gsutil_cp
#'
#' @param prior The name of gene sets used to annotate PCAGenomicSignatures. Currently
#' there are two available options.
#' \itemize{
#'     \item \code{C2} : MSigDB C2 (curated gene sets)
#'     \item \code{PLIERpriors} : bloodCellMarkersIRISDMAP, svmMarkers, and canonicalPathways
#' }
#' @param dir A path to the directory where PCAGenomicSginatures object will be downloaded.
#' @return PCAGenomicSignatures object in \code{.rds} format saved in the user
#' defined location.
#'
#' @export
getModel <- function(prior = c("C2", "PLIERpriors"), dir) {

  if (!prior %in% c("C2", "PLIERpriors")) {
    stop("Prior you entered isn't available yet.")
  }

  bucket <- "gs://pca_genomic_signatures"
  fname <- paste0("PCAmodel_", prior, ".rds")
  fpath <- file.path(bucket, fname)

  if (!dir.exists(dir)) {
    stop("Directory you entered doesn't exist.")
  } else {
    AnVIL::gsutil_cp(fpath, dir)
  }
}






# Download datasets used PCAGenomicSignatures manuscript
#
# getDataset <- function() {
# }
