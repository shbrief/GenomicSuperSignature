#' List the available RAVmodels
#' 
#' @param simplify Default is \code{TRUE}. If it is set to \code{FALSE}, the
#' additional metadata of different versions of RAVmodel
#' @return Under the default, this function will return a data frame with four
#' columns - prior, version, update, pkg_version. 
#' \itemize{
#'     \item \code{prior} : Different gene sets used for RAVmodel annotation.
#'     Currently, two are available - \code{C2} for MSigDB C2 (curated gene 
#'     sets), and \code{PLIERpriors} for bloodCellMarkersIRISDMAP, svmMarkers, 
#'     and canonicalPathways
#'     \item \code{version} : RAVmodel's version, which can be an input for 
#'     \code{version} argument of \code{\link{getModel}} function
#'     \item \code{update} : Date the RAVmodel is updated
#'     \item \code{pkg_version} : Compatible version of GenomicSuperSignature
#' }
#'
#' @examples 
#' availableRAVmodel()
#' 
#' @export
availableRAVmodel <- function(simplify = TRUE) {
    dir <- system.file("extdata", package = "GenomicSuperSignature")
    map <- utils::read.table(file.path(dir, "availableRAVmodel.csv"),
                             sep = ",", header = TRUE)

    if (isTRUE(simplify)) {
        available_ind <- which(map$gcp == "TRUE")
        map <- map[available_ind, 
                   c("prior", "version", "update", "pkg_version")]
    }
    
    return(map)
}
