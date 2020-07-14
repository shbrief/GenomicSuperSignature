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



### ----------------------------------------------
### Getter
### ----------------------------------------------

setGeneric("model", function(x, ...) standardGeneric("model"))
setMethod("model", "GenomicSignatures", function(x) {
    out <- assay(x)
    return(out)
})


setGeneric("geneSets", function(x, ...) standardGeneric("geneSets"))
setMethod("geneSets", "GenomicSignatures", function(x) {
    out <- S4Vectors::metadata(x)$geneSets
    return(out)
})


setGeneric("updateNote", function(x, ...) standardGeneric("updateNote"))
setMethod("updateNote", "GenomicSignatures", function(x) {
    out <- S4Vectors::metadata(x)$updateNote
    return(out)
})



### ----------------------------------------------
### Setter
### ----------------------------------------------

setGeneric("geneSets<-", function(x, ..., value) standardGeneric("geneSets<-"))
setMethod("geneSets<-", "GenomicSignatures", function(x, value) {
    S4Vectors::metadata(x)$geneSets <- value
    return(x)
})


setGeneric("updateNote<-", function(x, ..., value) standardGeneric("updateNote<-"))
setMethod("updateNote<-", "GenomicSignatures", function(x, value) {
    S4Vectors::metadata(x)$updateNote <- value
    return(x)
})
