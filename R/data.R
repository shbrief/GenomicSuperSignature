#' @name miniRAVmodel
#' @title RAVmodel from 536 studies, annotated with MSigDB C2
#' @docType data
#' @description
#' A object providing a miniature version of RAVmodel_C2 (PCAGenomicSignatures
#' object constructed from 536 studies and annotated with MSigDB C2).
#' @keywords data
#' @format PCAGenomicSignatures
#' @author Sehyun Oh \email{shbrief@gmail.com}
#'
"miniRAVmodel"


#' @name miniAllZ
#' @title Subset of allZ matrix constructed from 8 CRC training datasets
#' @docType data
#' @description
#' Eight colorectal cancer microarray datasets were used to build RAVmodel and
#' the intermediate file containing genes and top PCs from each dataset is named
#' as \code{allZ}. Hierarchical clustering result of \code{allZ} is saved as
#' \code{res_hcut}. For demonstration, we subset the \code{allZ} matrix with the
#' first 100 genes, which is named as \code{miniAllZ}.
#' @keywords data
#' @format A matrix with 100 genes and 160 PCs from 8 training datasets.
#' @author Sehyun Oh \email{shbrief@gmail.com}
#' @source https://github.com/shbrief/model_building/tree/main/RAVmodel_8CRC
#'
"miniAllZ"

#' @name res_hcut
#' @title Subset of allZ matrix constructed from 8 CRC training datasets
#' @docType data
#' @description
#' Eight colorectal cancer microarray datasets were used to build RAVmodel and
#' the intermediate file containing genes and top PCs from each dataset is named
#' as \code{allZ}. Hierarchical clustering result of \code{allZ} is saved as
#' \code{res_hcut}.
#' @keywords data
#' @format \code{hclust} object from \code{factoextra::hcut} function.
#' @author Sehyun Oh \email{shbrief@gmail.com}
#'
"res_hcut"

#' @name miniTCGA
#' @title Subset of TCGA-COAD and TCGA-BRCA RNA sequencing datasets
#' @docType data
#' @description
#' TCGA-COAD and TCGA-BRCA RNA sequencing data were acquired using
#' \code{GSEABenchmarkeR::loadEData} and log-transformed. Conversion from
#' EntrezID to gene symbol was done with \code{EnrichmentBrowser::idMap}. Only
#' 8 samples from each dataset are kept.
#' @keywords data
#' @format A list containing two SummarizedExperiment objects.
#' @author Sehyun Oh \email{shbrief@gmail.com}
#'
"miniTCGA"
