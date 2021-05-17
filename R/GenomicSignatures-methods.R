#' @include GenomicSignatures-class.R

#' @name GenomicSignatures-methods
#' @title Methods and accesors for \code{GenomicSignatures} object
#'
#' @description The default contents of \code{GenomicSignatures} object, with
#' a set of getter and setter generic functions, which extract either the
#' \code{assay}, \code{colData}, or \code{metadata} slots of a
#' \code{\link{GenomicSignatures-class}} object. When you create this object,
#' \code{colData$studies} should be populated before adding any information in
#' \code{trainingData} slot.
#'
#' @details
#' \itemize{
#'     \item assay(x) : RAVindex (= avgLoadings) containing genes x RAVs
#'     \item metadata(x) : Metadata associated with RAVindex building process
#'     \item colData(x) : Information on RAVs
#' }
#'
#' @section Setters:
#' Setter method values (i.e., \code{function(x) <- value}):
#' \itemize{
#'     \item metadata<- : Assign metadata
#'     \item coldata<- : Assign extra information associated with RAVs
#'     \item geneSets<- : A character vector containing the name of gene sets
#'     used to annotate average loadings
#'     \item updateNote<- : A character vector. Describes the main feature of a
#'     model construction
#' }
#'
#' @section Getters:
#' \itemize{
#'    \item RAVindex : Equivalent to \code{assays(x)$RAVindex}
#'    \item geneSets : Access the \code{metadata(x)$geneSets} slot
#'    \item updateNote : Access the \code{metadata(x)$updateNote} slot
#' }
#'
#' @param x A \code{GenomicSignatures} object
#' @param value See details.
#'
#' @return A GenomicSignatures object for the constructor
#'
#' @examples
#' data(miniRAVmodel)
#' miniRAVmodel
#'
#'
#' @aliases RAVindex colData metadata geneSets updateNote
#' geneSets<- updateNote<-
NULL


### ----------------------------------------------
### Getter
### ----------------------------------------------

#' @export
setGeneric("RAVindex", function(x) standardGeneric("RAVindex"))

#' @exportMethod RAVindex
#' @rdname GenomicSignatures-methods
setMethod("RAVindex", "GenomicSignatures", function(x) {
    out <- assay(x)
    return(out)
})

#' @export
setGeneric("geneSets", function(x) standardGeneric("geneSets"))

#' @exportMethod geneSets
#' @rdname GenomicSignatures-methods
setMethod("geneSets", "GenomicSignatures", function(x) {
    out <- metadata(x)$geneSets
    return(out)
})

#' @export
setGeneric("updateNote", function(x) standardGeneric("updateNote"))

#' @exportMethod updateNote
#' @rdname GenomicSignatures-methods
setMethod("updateNote", "GenomicSignatures", function(x) {
    out <- metadata(x)$updateNote
    return(out)
})



### ----------------------------------------------
### Setter
### ----------------------------------------------

#' @export
setGeneric("geneSets<-", function(x, value) standardGeneric("geneSets<-"))

#' @exportMethod geneSets<-
#' @rdname GenomicSignatures-methods
setMethod("geneSets<-", "GenomicSignatures", function(x, value) {
    metadata(x)$geneSets <- value
    return(x)
})

#' @export
setGeneric("updateNote<-", function(x, value) standardGeneric("updateNote<-"))

#' @exportMethod updateNote<-
#' @rdname GenomicSignatures-methods
setMethod("updateNote<-", "GenomicSignatures", function(x, value) {
    metadata(x)$updateNote <- value
    return(x)
})
