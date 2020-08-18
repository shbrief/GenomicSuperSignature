#' Validation result in heatmap format
#'
#' @param data An output matrix from \code{\link{validate}} function. Subset of
#' this matrix is plotted as a heatmap using \code{\link[ComplexHeatmap]{Heatmap}}
#' @param row_title A character string. Provide the row title.
#' @param num.out A number of highly validated PCclusters to output. Default is 5.
#' If any of the cutoff parameters are provided, \code{num.out} or the number of
#' filtered PCclusters, whichever smaller, will be chosen.
#' @param scoreCutoff A numeric value for the minimum correlation.
#' @param swCutoff A numeric value for the minimum average silhouette width.
#' @param clsizeCuoff A integer value for the minimum cluster size.
#' @param breaks A numeric vector of length 3. Number represents the values assigned
#' to three colors. Default is \code{c(0, 0.5, 1)}.
#' @param colors A character vector of length 3. Each represents the color assigned
#' to three breaks. Default is \code{c("white", "white smoke", "red")}.
#' @param whichPC An integer value between 1 and 8. PC number of your data to check
#' the validated signatures with. Under the default (\code{NULL}), it outputs top
#' scored signatures with any PC of your data.
#' @param ... any additional argument for \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return A heatmap with the validation result subset through the given conditions.
#'
#' @export
heatmapTable <- function(data, num.out = 5,
                         scoreCutoff = NULL, swCutoff = NULL, clsizeCutoff = NULL,
                         breaks = c(0, 0.5, 1),
                         colors = c("white", "white smoke", "red"),
                         column_title = NULL, row_title = NULL, whichPC = NULL, ...) {

  score_ind <- which(colnames(data) == "score")
  sw_ind <- which(colnames(data) == "sw")

  val <- validatedSignatures(data, num.out = num.out, whichPC = whichPC,
                             scoreCutoff = scoreCutoff, swCutoff = swCutoff, clsizeCutoff = clsizeCutoff)
  dat <- val[,score_ind,drop=FALSE] %>% t
  sw <- val[,sw_ind,drop=FALSE] %>% as.numeric

  # Define Heatmap Table Size
  if (is.null(column_title)) {
    column_title <- character(0)
    hmHeight <- nrow(dat) + 2.5
  } else {
    column_title <- column_title
    hmHeight <- nrow(dat) + 4
  }

  if (is.null(row_title)) {
    row_title <- character(0)
    hmWidth <- ncol(dat) + 2
  } else {
    row_title <- row_title
    hmWidth <- ncol(dat) + 3
  }

  # Heatmap annotation representing Silhouette width of selected PCclusters
  ha <- ComplexHeatmap::HeatmapAnnotation(avg.sw = ComplexHeatmap::anno_barplot(sw, baseline = 0,
                                                                                gp = grid::gpar(fill = ifelse(sw > 0, "red", "blue"), alpha = 0.5)),
                                          annotation_name_side = "left",
                                          annotation_name_gp = grid::gpar(fontsize = 10))

  # Draw Heatmap
  ComplexHeatmap::Heatmap(dat, col = circlize::colorRamp2(breaks, colors), name = "Corr",
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid::grid.text(sprintf("%.2f", dat[i, j]), x, y, gp = grid::gpar(fontsize = 8))
                          },
                          row_title = row_title,
                          row_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
                          row_names_gp = grid::gpar(fontsize = 10),
                          row_names_side = "left",
                          row_names_max_width = unit(0.5, "cm"),
                          column_title = column_title,
                          column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
                          column_names_side = "bottom",
                          column_names_gp = grid::gpar(fontsize = 10),
                          column_names_rot = 90,
                          column_names_max_height = unit(1, "cm"),
                          heatmap_width = unit(hmWidth, "cm"),
                          heatmap_height = unit(hmHeight, "cm"),
                          top_annotation = ha, ...
  )}
