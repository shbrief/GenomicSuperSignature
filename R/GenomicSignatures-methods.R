#' @include PCAGenomicSignatures-methods.R


### ----------------------------------------------
### Getter
### ----------------------------------------------

#' @export
setGeneric("model", function(x) standardGeneric("model"))

#' @exportMethod model
#' @rdname PCAGenomicSignatures-methods
setMethod("model", "GenomicSignatures", function(x) {
    out <- assay(x)
    return(out)
})

#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))

#' @exportMethod metadata
#' @rdname PCAGenomicSignatures-methods
setMethod("metadata", "GenomicSignatures", function(x) {
    out <- S4Vectors::metadata(x)
    return(out)
})

#' @export
setGeneric("geneSets", function(x) standardGeneric("geneSets"))

#' @exportMethod geneSets
#' @rdname PCAGenomicSignatures-methods
setMethod("geneSets", "GenomicSignatures", function(x) {
    out <- S4Vectors::metadata(x)$geneSets
    return(out)
})


#' @export
setGeneric("updateNote", function(x) standardGeneric("updateNote"))

#' @exportMethod updateNote
#' @rdname PCAGenomicSignatures-methods
setMethod("updateNote", "GenomicSignatures", function(x) {
    out <- S4Vectors::metadata(x)$updateNote
    return(out)
})



### ----------------------------------------------
### Setter
### ----------------------------------------------

#' @export
setGeneric("geneSets<-", function(x, value) standardGeneric("geneSets<-"))

#' @exportMethod geneSets<-
#' @rdname PCAGenomicSignatures-methods
setMethod("geneSets<-", "GenomicSignatures", function(x, value) {
    S4Vectors::metadata(x)$geneSets <- value
    return(x)
})


#' @export
setGeneric("updateNote<-", function(x, value) standardGeneric("updateNote<-"))

#' @exportMethod updateNote<-
#' @rdname PCAGenomicSignatures-methods
setMethod("updateNote<-", "GenomicSignatures", function(x, value) {
    S4Vectors::metadata(x)$updateNote <- value
    return(x)
})
