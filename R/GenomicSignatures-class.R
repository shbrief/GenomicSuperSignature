### ==============================================
### Class
### ==============================================
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @exportClass GenomicSignatures
GenomicSignatures <- setClass("GenomicSignatures",
                              contains = c("SummarizedExperiment", "VIRTUAL")
)

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
