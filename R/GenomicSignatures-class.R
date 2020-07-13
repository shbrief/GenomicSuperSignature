### ==============================================
### GenomicSignatures-Class
### ==============================================
#' GenomicSignatures
#' @title Construct \code{GenomicSignatures} object
#' @name GenomicSignatures-class
#' @import SummarizedExperiment
#' @return A \code{GenomicSignatures} object
#' @exportClass GenomicSignatures
GenomicSignatures <- setClass("GenomicSignatures",
                              contains = c("SummarizedExperiment", "VIRTUAL")
)

#' @name GenomicSignatures-methods
#' @title Accessing and modifying information in GenomicSignatures
#'
#' @description A set of accessor and setter generic functions to extract
#' either the \code{assay}, \code{colData}, or \code{metadata} slots of a
#' \code{\link{GenomicSignatures}} object
#'
#' @param x A \code{GenomicSignatures} object
#' @param value See details.
NULL


### ----------------------------------------------
### Getter
### ----------------------------------------------
#' @export
setGeneric("model", function(x, ...) standardGeneric("model"))

#' @exportMethod model
#' @rdname GenomicSignatures-methods
setMethod("model", "GenomicSignatures", function(x) {
    out <- assay(x)
    return(out)
})

#' @export
setGeneric("geneSets", function(x, ...) standardGeneric("geneSets"))

#' @exportMethod geneSets
#' @rdname GenomicSignatures-methods
setMethod("geneSets", "GenomicSignatures", function(x) {
    out <- metadata(x)$geneSets
    return(out)
})

#' @export
setGeneric("updateNote", function(x, ...) standardGeneric("updateNote"))

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
setGeneric("geneSets<-", function(x, ..., value) standardGeneric("geneSets<-"))

#' @importMethodsFrom S4Vectors metadata
#' @exportMethod geneSets<-
#' @rdname GenomicSignatures-methods
setMethod("geneSets<-", "GenomicSignatures", function(x, value) {
    metadata(x)$geneSets <- value
    return(x)
})

#' @export
setGeneric("updateNote<-", function(x, ..., value) standardGeneric("updateNote<-"))

#' @exportMethod updateNote<-
#' @rdname GenomicSignatures-methods
setMethod("updateNote<-", "GenomicSignatures", function(x, value) {
    metadata(x)$updateNote <- value
    return(x)
})
