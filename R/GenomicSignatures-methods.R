#' @include GenomicSignatures-class.R

#' @name GenomicSignatures-methods
#' @title Methods and accesors for \code{GenomicSignatures} object
#'
#' @description The default contents of \code{GenomicSignatures} object, with
#' a set of accessor and setter generic functions, which extract either the \code{assay},
#' \code{colData}, or \code{metadata} slots of a \code{\link{GenomicSignatures-class}}
#' object. When you create this object, \code{colData$studies} should be populated
#' before adding any information in \code{trainingData} slot
#'
#' @details
#' \itemize{
#'     \item assay(x) : PCAmodel (= avgLoading) containing genes x PCclusters
#'     \item metadata(x) : Metadata associated with the PCAmodel buildling process
#'     \item colData(x) : Information on PCclusters
#' }
#'
#' @section Setters:
#' Setter method values (i.e., \code{function(x) <- value}):
#' \itemize{
#'     \item geneSets<- : A character vector containing the name of gene sets used
#'     to annotate average loadings
#'     \item updateNote<- : A character vetor. Describes the main feature of a model construction
#' }
#'
#' @section Accessors:
#' \itemize{
#'    \item model : Equivalent to the \code{assays(x)$model} accessor for convenience
#'    \item geneSets : Access the \code{metadata(x)$geneSets} slot
#'    \item updateNote : Access the \code{metadata(x)$updateNote} slot
#' }
#'
#' @param x A \code{GenomicSignatures} object
#' @param value See details.
#'
#' @aliases model colData metadata geneSets updateNote geneSets<- updateNote<-
NULL


### ----------------------------------------------
### Getter
### ----------------------------------------------

#' @export
setGeneric("model", function(x) standardGeneric("model"))

#' @exportMethod model
#' @rdname GenomicSignatures-methods
setMethod("model", "GenomicSignatures", function(x) {
    out <- assay(x)
    return(out)
})

#' @export
setGeneric("colData", function(x) standardGeneric("colData"))

#' @exportMethod colData
#' @rdname GenomicSignatures-methods
setMethod("colData", "GenomicSignatures", function(x) {
    getElement(x, "colData")
})

#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))

#' @exportMethod metadata
#' @rdname GenomicSignatures-methods
setMethod("metadata", "GenomicSignatures", function(x)
    getElement(x, "metadata"))

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

#' @exportMethod metadata<-
#' @rdname GenomicSignatures-methods
setReplaceMethod("metadata", c("GenomicSignatures", "ANY"),
                 function(x, ..., value) {
                     slot(x, "metadata") <- value
                     return(x)
                 })


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
