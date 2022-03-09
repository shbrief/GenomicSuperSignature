s2p_get_cache <- function(cache = tools::R_user_dir("GenomicSuperSignature",
                                                    "cache")) {
  BiocFileCache::BiocFileCache(cache=cache)
}

#' @importFrom BiocFileCache bfcneedsupdate bfcdownload bfcadd bfcquery bfcrpath
#'
s2p_cached_url <- function(url, rname = url, ask_on_update=FALSE,
                           ...) {
  bfc = s2p_get_cache()
  bfcres = bfcquery(bfc,rname,'rname')

  rid = bfcres$rid
  # Not found
  fileage = 0
  if(!length(rid)) {
    rid = names(bfcadd(bfc, rname, url))
  }
  # if needs update, do the download
  if(bfcneedsupdate(bfc, rid)) {
    bfcdownload(bfc, rid, ask=FALSE, ...)
    print("downloading")
  }
  bfcrpath(bfc, rids = rid)
}




#' Download a PCAGenomicSignatures model
#'
#' @param prior The name of gene sets used to annotate PCAGenomicSignatures.
#' Currently  there are two available options.
#' \itemize{
#'     \item \code{C2} : MSigDB C2 (curated gene sets)
#'     \item \code{PLIERpriors} : bloodCellMarkersIRISDMAP, svmMarkers, and
#'     canonicalPathways
#' }
#' @param load Default is \code{TRUE}. If it's set to \code{FALSE}, the models
#' are just downloaded to cache but not loaded into memory.
#' @param version Default is \code{latest}. Available versions are listed in
#' \code{version} column of \code{\link{availableRAVmodel()}} output.
#'
#' @return File cache location or PCAGenomicSignatures object loaded from it.
#'
#' @examples
#' z = getModel("C2")
#'
#' @export
getModel <- function(prior = c("C2", "PLIERpriors"), 
                     version = "latest", load = TRUE) {

    dir <- system.file("extdata", package = "GenomicSuperSignature")
    models <- read.table(file.path(dir, "availableRAVmodel.csv"), 
                         sep = ",", header=TRUE)
    models <- models[which(models$gcp == TRUE),] ## Only RAVmodels on GCP
    
    ## Check the input validity
    if (!prior %in% unique(models$prior)) {
      stop("Prior you entered isn't available.")
    }
    
    prior_ind <- which(models$prior == prior)
    if (!version %in% unique(models$version[prior_ind])) {
      stop("Version you entered doesn't exist.")
    }
    
    request_ind <- which(models$prior == prior & models$version == version)
    request_model <- paste0(models$fname[request_ind], ".rds")

    bucket_name <- "genomic_super_signature"
    fpath <- file.path("https://storage.googleapis.com",
                       bucket_name, request_model)

    fpath <- s2p_cached_url(fpath)
    if (isTRUE(load)) {
      model <- readRDS(fpath)
      return(model)
    } else {return(fpath)}
}



# Download datasets used PCAGenomicSignatures manuscript
#
# getDataset <- function() {
# }
