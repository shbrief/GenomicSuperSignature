% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{miniTCGA}
\alias{miniTCGA}
\title{Subset of TCGA-COAD and TCGA-BRCA RNA sequencing datasets}
\format{
A list containing two SummarizedExperiment objects.
}
\usage{
miniTCGA
}
\description{
TCGA-COAD and TCGA-BRCA RNA sequencing data were acquired using
\code{GSEABenchmarkeR::loadEData} and log-transformed. Conversion from
EntrezID to gene symbol was done with \code{EnrichmentBrowser::idMap}. Only
8 samples from each dataset are kept.
}
\author{
Sehyun Oh \email{shbrief@gmail.com}
}
\keyword{data}
