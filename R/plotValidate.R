#' Plot validation results in an interactive graph
#'
#' There are three main information on the graph:
#' \itemize{
#'     \item x-axis : Pearson correlation coefficient. Higher value means that
#'     test dataset and RAV is more tightly associated with.
#'     \item y-axis : Silhouette width representing the quality of RAVs.
#'     \item size : The number of studies in each RAV. (= cluster size)
#'     \item color : Test dataset's PC number that validate each RAV. Because we
#'     used top 8 PCs of the test dataset, there are 8 categories.
#' }
#'
#' @import ggplot2
#'
#' @param val_all Output from \link{validate} function.
#' @param swFilter If \code{swFilter=TRUE}, only RAV above the cutoff, defined
#' through \code{minSilhouetteWidth} argument will be plotted. Default is
#' \code{swFilter=FALSE}
#' @param minClusterSize The minimum size of clusters to be included in the
#' plotting. Default value is 2, so any single-element clusters are excluded.
#' @param minSilhouetteWidth A minimum average silhouette width to be plotted.
#' Only effective under \code{swFilter=TRUE} condition. Default is 0.
#' @param interactive If set to \code{TRUE}, the output will be interactive
#' plot. Default is \code{FALSE}.
#' @param minClSize The minimum number of PCs in the clusters you want.
#' @param maxClSize The maximum number of PCs in the clusters you want.
#' @param colorPalette Default is \code{Dark2}. For other color options, please
#' check \code{\link[ggplot2]{scale_color_brewer}}
#'
#' @return a ggplot object
#'
#' @examples
#' data(miniRAVmodel)
#' library(bcellViper)
#' data(bcellViper)
#' val_all <- validate(dset, miniRAVmodel)
#' plotValidate(val_all)
#'
#' @export
plotValidate <- function(val_all, minClusterSize = 2, swFilter = FALSE,
                         minSilhouetteWidth = 0, interactive = FALSE,
                         minClSize = NULL, maxClSize = NULL,
                         colorPalette = "Dark2") {

    # Validate inputs
    stopifnot(length(swFilter) == 1L, !is.na(swFilter), is.logical(swFilter))
    stopifnot(length(interactive) == 1L,
              !is.na(interactive), is.logical(interactive))

    # If the validation result contains all PCs (val_all from `level = "all"`)
    if (identical(colnames(val_all), paste0("PC", seq_len(8)))) {
        stop("'val_all' input should be created by `validate` function with
         `level = \"max\"`, not `level = \"all\"`.")
    }

    # If the validation result is from the list of datasets
    score_ind <- which(colnames(val_all) == "score")
    if (length(score_ind) == 0) {
        stop("'val_all' input should be from one dataset,
             not from multiple datasets")
    }

    ## Binding the variables from val_all locally to the function
    cl_size <- sw <- score <- cl_num <- PC <- NULL

    val_all <-  as.data.frame(val_all)
    clSizeFiltered <- dplyr::filter(val_all, cl_size >= minClusterSize)

    if (swFilter) {
        filtered <- dplyr::filter(clSizeFiltered, sw >= minSilhouetteWidth)
    } else {filtered <- clSizeFiltered}

    ind <- which(colnames(filtered) == "PC")
    filtered[,ind] <- factor(filtered[,ind])

    p <- ggplot(filtered, aes(sw, score,
                              size = cl_size,
                              label = cl_num,
                              color = PC)) +
        geom_point(alpha = 0.5) +
        scale_size(range = c(0.5, 8), name = "Cluster Size") +
        scale_color_brewer(palette = colorPalette) +
        theme_bw() +
        xlab("Average silhouette width of each cluster") +
        ylab("Validation score")

    if (interactive) {plotly::ggplotly(p)} else {p}
}

