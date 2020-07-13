### ==============================================
### PCAGenomicSignatures-Class
### ==============================================
#' PCAGenomicSignatures
#' @import methods SummarizedExperiment S4Vectors
#' @include GenomicSignatures-class.R
#' @exportClass PCAGenomicSignatures
#' @slot trainingData Containing metadata of each studies used in model building
.PCAGenomicSignatures <- setClass("PCAGenomicSignatures",
                                  slots = representation(
                                      trainingData = "DataFrame"
                                  ),
                                  contains = "GenomicSignatures"
)

#' Formatting PCcluster name
#'
#' Keep the name with 'k + cluster number + number of PCs + number of unique studies'
#' info during the model construction to make it easy to keep track of them, but at the
#' PCAGenomicSignatures-class object building step, covert them into a short version,
#' mPC (mean of PCs).
#'
#' @param x PCAGenomicSignature-object
#' @param ... Additional arguments for supporting functions.
.PCclusterName <- function(x, ...) {
    x@colData$PCcluster <- rownames(x@colData)
    rownames(x@colData) <- paste0("PCcluster", 1:ncol(x))
    return(x)
}


### ==============================================
### PCAGenomicSignatures Constructor
### ==============================================
#' @name PCAGenomicSignatures
#' @title Construct \code{PCAGenomicSignatures} object
#' @description The default contents of \code{PCAGenomicSignatures} object, with a set
#' of accessor and setter generic functions. When you create this object, \code{colData(x)$studies}
#' should be populated before adding any information in \code{trainingData(x)} slot
#' \itemize{
#'     \item assay(x) : PCAmodel (= avgLoading) containing genes x PCclusters
#'     \item metadata(x)$cluster : A vector of integers (from 1:k) indicating the
#'     cluster to which each point is allocated. Check \code{stats::kmeans} function.
#'     \item metadata(x)$size : The number of points in each cluster. Check \code{stats::kmeans} function.
#'     \item metadata(x)$iter : The number of (outer) iterations. Check \code{stats::kmeans} function.
#'     \item metadata(x)$k : The number of PCclusters.
#'     \item metadata(x)$n : The number of top PCs from each dataset.
#'     \item metadata(x)$geneSets : Name of the prior gene sets used to annotate average loadings.
#'     \item colData(x)$studies : A list of character vectors containing studies contributing to each PC cluster.
#'     \item colData(x)$gsea : A list of \code{gseaResult} objects. Build using \code{clusterProfiler::GSEA} function.
#' }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @param ... Additional arguments for supporting functions
#' @export
PCAGenomicSignatures <- function(...)
{
    se <- SummarizedExperiment(...)
    gs <- .PCAGenomicSignatures(se)
    .PCclusterName(gs)
}


#' @name PCAGenomicSignatures-methods
#' @title Accessing and modifying information in PCAGenomicSignatures
#'
#' @description A set of accessor and setter generic functions to extract
#' either the \code{assay}, \code{colData}, or \code{metadata} slots of a
#' \code{\link{PCAGenomicSignatures}} object
#'
#' @section Setters:
#' Setter method values (i.e., \code{function(x) <- value}):
#' \itemize{
#'     \item geneSets<- : A character vector containing the name of gene sets used
#'     to annotate average loadings
#'     \item studies<- : A list of character vectors containing gene sets used to annotate average loadings
#'     \item gsea<- : A list of \code{gseaResult} objects.
#'     \item metadata<- : A \code{list} object of metadata
#'     \item `$<-` : A vector to replace the indicated column in \code{colData}
#' }
#'
#' @section Accessors:
#' All the accessors inherited from \code{SummarizedExperiment} are available and
#' the additional accessors for \code{PCAGenomicSignatures} specific data are listed
#' below.
#' \itemize{
#'    \item model : Equivalent to the \code{assays(x)$model} accessor for convenience
#'    \item geneSets : Access the \code{metadata(x)$geneSets} slot
#'    \item studies : Access the \code{colData(x)$studies} slot
#'    \item gsea : Access the \code{colData(x)$gsea}
#'    \item `$` : Access a column in \code{colData}
#'    \item trainingData : Access the \code{trainingData} slot
#'    \item mesh : Access the \code{trainingData(x)$MeSH} slot
#'    \item PCAsummary : Access the \code{trainingData(x)$PCAsummary} slot
#' }
#'
#' @param x A \code{PCAGenomicSignatures} object
#' @param value See details.
NULL


### ==============================================
### Setter
### ==============================================

### PCAGenomicSignatures
#' @export
setGeneric("studies<-", function(x, ..., value) standardGeneric("studies<-"))

#' @import dplyr
#' @exportMethod studies<-
#' @rdname PCAGenomicSignatures-methods Studies associated with each cluster
setMethod("studies<-", "PCAGenomicSignatures", function(x, value) {
    x@colData$studies <- value
    allStudies <- unlist(value) %>% unique
    x@trainingData <- DataFrame(row.names = allStudies)

    # validObject(x)
    return(x)
})

#' @export
setGeneric("silhouetteWidth<-", function(x, ..., value) standardGeneric("silhouetteWidth<-"))

#' @import dplyr
#' @exportMethod silhouetteWidth<-
#' @rdname PCAGenomicSignatures-methods Average Silhouette width of each cluster
setMethod("silhouetteWidth<-", "PCAGenomicSignatures", function(x, value) {
    x@colData$silhouetteWidth <- value
    return(x)
})


#' @export
setGeneric("gsea<-", function(x, ..., value) standardGeneric("gsea<-"))

#' @exportMethod gsea<-
#' @rdname PCAGenomicSignatures-methods GSEA output of each cluster
setMethod("gsea<-", "PCAGenomicSignatures", function(x, value) {

    allValue <- vector(mode = "list", length = length(colnames(x)))
    names(allValue) <- colnames(x)
    ind = which(x@colData$PCcluster %in% names(value))
    allValue[ind] = value

    x@colData$gsea = allValue
    return(x)
})

#' @export
setGeneric("trainingData<-", function(x, ..., value) standardGeneric("trainingData<-"))

#' @import dplyr
#' @exportMethod trainingData<-
#' @rdname PCAGenomicSignatures-methods Metadata, such as PCAsummary and MeSH, associated
#' with each study used for model
setMethod("trainingData<-", "PCAGenomicSignatures", function(x, value) {
    x@trainingData <- value
    # validObject(x)
    return(x)
})

#' @export
setGeneric("mesh<-", function(x, ..., value) standardGeneric("mesh<-"))

#' @exportMethod mesh<-
#' @rdname PCAGenomicSignatures-methods MeSH terms assigned to each study
setMethod("mesh<-", "PCAGenomicSignatures", function(x, value) {
    trainingData(x)$MeSH = NA
    for (i in seq_along(value)) {
        ind = which(rownames(trainingData(x)) == names(value[i]))
        trainingData(x)$MeSH[ind] = value[i]
    }
    names(trainingData(x)$MeSH) = rownames(trainingData(x))
    return(x)
})


#' @export
setGeneric("PCAsummary<-", function(x, ..., value) standardGeneric("PCAsummary<-"))

#' @exportMethod PCAsummary<-
#' @rdname PCAGenomicSignatures-methods PCA output summary including variance explained
setMethod("PCAsummary<-", "PCAGenomicSignatures", function(x, value) {
    trainingData(x)$PCAsummary = NA
    for (i in seq_along(value)) {
        ind = which(rownames(trainingData(x)) == names(value[i]))
        trainingData(x)$PCAsummary[ind] = value[i]
    }
    names(trainingData(x)$PCAsummary) = rownames(trainingData(x))
    return(x)
})








### ==============================================
### Getter
### ==============================================

### PCAGenomicSignatures
#' @export
setGeneric("studies", function(x, ...) standardGeneric("studies"))

#' @exportMethod studies
#' @rdname PCAGenomicSignatures-methods Access studies in each cluster
setMethod("studies", "PCAGenomicSignatures", function(x) {
    out <- x@colData$studies
    return(out)
})


#' @export
setGeneric("silhouetteWidth", function(x, ...) standardGeneric("silhouetteWidth"))

#' @exportMethod silhouetteWidth
#' @rdname PCAGenomicSignatures-methods Access average silhouette width of each cluster
setMethod("silhouetteWidth", "PCAGenomicSignatures", function(x) {
    out <- x@colData$silhouetteWidth
    return(out)
})


#' @export
setGeneric("gsea", function(x, ...) standardGeneric("gsea"))

#' @exportMethod gsea
#' @rdname PCAGenomicSignatures-methods Acess GSEA output of each cluster
setMethod("gsea", "PCAGenomicSignatures", function(x) {
    out <- x@colData$gsea
    return(out)
})

#' @export
setGeneric("trainingData", function(x, ...) standardGeneric("trainingData"))

#' @exportMethod trainingData
#' @rdname PCAGenomicSignatures-methods Access metadata, such as PCAsummary and MeSH, of clusters
setMethod("trainingData", "PCAGenomicSignatures", function(x) {
    out <- x@trainingData
    return(out)
})

#' @export
setGeneric("mesh", function(x, ...) standardGeneric("mesh"))

#' @exportMethod mesh
#' @rdname PCAGenomicSignatures-methods Access MeSH terms assigned to each study
setMethod("mesh", "PCAGenomicSignatures", function(x) {
    out <- x@trainingData$MeSH
    return(out)
})

#' @export
setGeneric("PCAsummary", function(x, ...) standardGeneric("PCAsummary"))

#' @exportMethod PCAsummary
#' @rdname PCAGenomicSignatures-methods Access PCA summary of each study
setMethod("PCAsummary", "PCAGenomicSignatures", function(x) {
    out <- x@trainingData$PCAsummary
    return(out)
})


### ==============================================
### Show method
### ==============================================
#' @import dplyr
#' @import S4Vectors
#' @export
setMethod("show", "PCAGenomicSignatures", function(object) {
    callNextMethod()

    colnames <- colnames(trainingData(object))
    coolcat("trainingData(%d): %s\n", colnames)

    rownames <- rownames(trainingData(object))
    coolcat("trainingData names(%d): %s\n", rownames)
})

