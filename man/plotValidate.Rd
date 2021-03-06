% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotValidate.R
\name{plotValidate}
\alias{plotValidate}
\title{Plot validation results in an interactive graph}
\usage{
plotValidate(
  val_all,
  minClusterSize = 2,
  swFilter = FALSE,
  minSilhouetteWidth = 0,
  interactive = FALSE,
  minClSize = NULL,
  maxClSize = NULL,
  colorPalette = "Dark2"
)
}
\arguments{
\item{val_all}{Output from \link{validate} function.}

\item{minClusterSize}{The minimum size of clusters to be included in the
plotting. Default value is 2, so any single-element clusters are excluded.}

\item{swFilter}{If \code{swFilter=TRUE}, only RAV above the cutoff, defined
through \code{minSilhouetteWidth} argument will be plotted. Default is
\code{swFilter=FALSE}}

\item{minSilhouetteWidth}{A minimum average silhouette width to be plotted.
Only effective under \code{swFilter=TRUE} condition. Default is 0.}

\item{interactive}{If set to \code{TRUE}, the output will be interactive
plot. Default is \code{FALSE}.}

\item{minClSize}{The minimum number of PCs in the clusters you want.}

\item{maxClSize}{The maximum number of PCs in the clusters you want.}

\item{colorPalette}{Default is \code{Dark2}. For other color options, please
check \code{\link[ggplot2]{scale_color_brewer}}}
}
\value{
a ggplot object
}
\description{
There are three main information on the graph:
\itemize{
    \item x-axis : Pearson correlation coefficient. Higher value means that
    test dataset and RAV is more tightly associated with.
    \item y-axis : Silhouette width representing the quality of RAVs.
    \item size : The number of studies in each RAV. (= cluster size)
    \item color : Test dataset's PC number that validate each RAV. Because we
    used top 8 PCs of the test dataset, there are 8 categories.
}
}
\examples{
data(miniRAVmodel)
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, miniRAVmodel)
plotValidate(val_all)

}
