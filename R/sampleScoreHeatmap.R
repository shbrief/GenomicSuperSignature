#' Plot heatmap of the sample scores
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar
#'
#' @param score An output from \code{\link{calculateScore}} function, which is a matrix
#' with samples (row) and PrcompClusters (column) If it is a simple vector, it
#' will be converted to a one-column matrix.
#' @param dataName Title on the row. The name of the dataset to be scored.
#' @param modelName Title on the column. The PCAmodel used for scoring.
#' @param cluster_rows A logical. Under the default (\code{TRUE}), rows will be clustered.
#' @param cluster_columns A logical. Under the default (\code{TRUE}), columns will be clustered.
#' @param show_row_names Whether show row names. Default is \code{TRUE}, showing the row name.
#' @param show_column_names Whether show column names. Default is \code{TRUE}, showing the column name.
#' @param row_names_gp Graphic parameters for row names. The default is 0.7.
#' @param column_names_gp Graphic parameters for column names. The default is 5.
#' @param ... Any additional argument for \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return A heatmap of the sample score. Rows represent samples and columns
#' represent PCclusters.
#'
#' @examples
#' data(miniPCAmodel)
#' library(bcellViper)
#' data(bcellViper)
#' score <- calculateScore(dset, miniPCAmodel)
#' sampleScoreHeatmap(score, dataName = "bcellViper", modelName = "miniPCAmodel")
#'
#' @export
sampleScoreHeatmap <- function(score, dataName, modelName,
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               show_row_names = TRUE,
                               show_column_names = TRUE,
                               row_names_gp = 0.7,
                               column_names_gp = 5, ...) {
    ComplexHeatmap::Heatmap(score,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            row_title = dataName,
                            row_title_gp = gpar(fontsize = 11, fontface = "bold"),
                            column_title = modelName,
                            column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                            row_names_gp = gpar(fontsize = row_names_gp),
                            column_names_gp = gpar(fontsize = column_names_gp),
                            show_row_names = show_row_names,
                            show_column_names = show_column_names,
                            heatmap_width = unit(10, "cm"),
                            heatmap_legend_param = list(title = "score"), ...)
}
