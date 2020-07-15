### ==============================================
### GenomicSignatures-Class
### ==============================================
#' GenomicSignatures
#' @name GenomicSignatures-class
#' @title Construct GenomicSignatures-class object
#'
#' @description GenomicSignatures is a virtual class inherited from SummarizedExperiment
#' and hosts GenomicSignatures models built from different dimensional reduction methods.
#' Currently, PCA-based model, called PCAGenomicSignatures, is available.
#'
#' @param x A \code{GenomicSignatures-class} object
#' @param value See details.
#'
#' @import SummarizedExperiment
#' @docType class
#' @exportClass GenomicSignatures
setClass("GenomicSignatures",
         contains = c("SummarizedExperiment", "VIRTUAL")
)

### ==============================================
### PCAGenomicSignatures-Class
### ==============================================
#' PCAGenomicSignatures
#' @name PCAGenomicSignatures-class
#' @title Construct PCAGenomicSignatures-class object
#'
#' @description PCA-based GenomicSignatures.
#'
#' @slot trainingData A \code{\link[S4Vectors]{DataFrame}} class object for metadata
#' associated with training data
#' @param x A \code{GenomicSignatures-class} object
#' @param value See details.
#'
#' @docType class
#' @exportClass PCAGenomicSignatures
.PCAGenomicSignatures <- setClass("PCAGenomicSignatures",
         slots = representation(
             trainingData = "DataFrame"
         ),
         contains = "GenomicSignatures"
)
