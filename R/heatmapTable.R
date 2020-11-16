.heatmapTableMultiStudies <- function(val_all, scoreCutoff = scoreCutoff,
                                      breaks = breaks, colors = colors,
                                      column_title = column_title,
                                      row_title = row_title, ...) {

  dat <- validatedSignatures(val_all, scoreCutoff = scoreCutoff)

  # Define Heatmap Table Size
  if (is.null(column_title)) {
    column_title <- character(0)
    hmHeight <- nrow(dat) + 1
  } else {
    column_title <- column_title
    hmHeight <- nrow(dat) + 2
  }

  if (is.null(row_title)) {
    row_title <- character(0)
    hmWidth <- ncol(dat) + 1
  } else {
    row_title <- row_title
    hmWidth <- ncol(dat) + 2
  }

  # Draw Heatmap
  ComplexHeatmap::Heatmap(dat, col = circlize::colorRamp2(breaks, colors), name = "Corr",
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid::grid.text(sprintf("%.2f", dat[i, j]), x, y, gp = grid::gpar(fontsize = 8))
                          },
                          show_heatmap_legend = FALSE,
                          row_title = row_title,
                          row_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
                          row_names_gp = grid::gpar(fontsize = 10),
                          row_names_side = "right",
                          row_names_max_width = unit(0.5, "cm"),
                          column_title = column_title,
                          column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
                          column_names_side = "bottom",
                          column_names_gp = grid::gpar(fontsize = 10),
                          column_names_rot = 90,
                          column_names_max_height = unit(1, "cm"),
                          heatmap_width = unit(hmWidth, "cm"),
                          heatmap_height = unit(hmHeight, "cm"), ...
  )
}

#' Validation result in heatmap format
#'
#' This function subsets \code{\link{validate}} outputs with different criteria
#' and visualize it in a heatmap-like table.
#'
#' @param val_all An output matrix from \code{\link{validate}} function with the
#' parameter \code{level = "max"}. Subset of this matrix is plotted as a heatmap
#' using \code{\link[ComplexHeatmap]{Heatmap}}
#' @param ind An integer vector. If this parameter is provided, the other parameters,
#' \code{num.out, scoreCutoff, swCutoff, clsizeCutoff} will be ignored and the heatmap
#' table containing only the provided index will be printed.
#' @param num.out A number of highly validated RAVs to output. Default is 5.
#' If any of the cutoff parameters are provided, \code{num.out} or the number of
#' filtered RAVs, whichever smaller, will be chosen.
#' @param scoreCutoff A numeric value for the minimum correlation. If \code{val_all}
#' input is from multiple studies, the default is 0.7 and this is the only cutoff
#' criteria considred: \code{swCutoff} and \code{clsizeCutoff} will be ignored.
#' @param swCutoff A numeric value for the minimum average silhouette width.
#' @param clsizeCutoff A integer value for the minimum cluster size.
#' @param breaks A numeric vector of length 3. Number represents the values assigned
#' to three colors. Default is \code{c(0, 0.5, 1)}.
#' @param colors A character vector of length 3. Each represents the color assigned
#' to three breaks. Default is \code{c("white", "white smoke", "red")}.
#' @param column_title A character string. Provide the column title.
#' @param row_title A character string. Provide the row title.
#' @param whichPC An integer value between 1 and 8. PC number of your data to check
#' the validated signatures with. Under the default (\code{NULL}), it outputs top
#' scored signatures with any PC of your data.
#' @param ... any additional argument for \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return A heatmap displaying the subset of the validation result that met the
#' given cutoff criteria. If \code{val_all} input is from a single dataset, the
#' output heatmap will contain both score and average silhouette width for each
#' cluster.
#'
#' If \code{val_all} input is from multiple studies, the output heatmap's rows will
#' represent each study and the columns will be RAVs, which meet \code{scoreCutoff}
#' for any of the input studies.
#'
#' @examples
#' data(miniPCAmodel)
#' library(bcellViper)
#' data(bcellViper)
#' val_all <- validate(dset, miniPCAmodel)
#' heatmapTable(val_all, swCutoff = 0)
#'
#' @export
heatmapTable <- function(val_all, ind = NULL, num.out = 5,
                         scoreCutoff = NULL, swCutoff = NULL, clsizeCutoff = NULL,
                         breaks = c(0, 0.5, 1),
                         colors = c("white", "white smoke", "red"),
                         column_title = NULL, row_title = NULL, whichPC = NULL, ...) {

  score_ind <- which(colnames(val_all) == "score")
  sw_ind <- which(colnames(val_all) == "sw")

  # If the validation result contains all PCs (`validate` with `level = "all"`)
  if (identical(colnames(val_all), paste0("PC", 1:8))) {
    stop("'val_all' input should be created by `validate` function with
         `level = \"max\"`, not `level = \"all\"`.")
  }

  # If the validation result is from the list of datsets
  if (length(score_ind) == 0) {
    res <- .heatmapTableMultiStudies(val_all, scoreCutoff = scoreCutoff,
                                     breaks = breaks, colors = colors,
                                     column_title = column_title,
                                     row_title = row_title)
    return(res)
    stop
  }

  if (is.null(ind)) {
    val <- validatedSignatures(val_all, num.out = num.out, whichPC = whichPC, indexOnly = FALSE,
                               scoreCutoff = scoreCutoff, swCutoff = swCutoff, clsizeCutoff = clsizeCutoff)
    dat <- val[,score_ind,drop=FALSE] %>% t
    sw <- val[,sw_ind,drop=FALSE] %>% as.numeric
  } else {
    row_ind <- which(rownames(val_all) %in% paste0("RAV", ind))
    val <- val_all[row_ind,]
    dat <- val[,score_ind,drop=FALSE] %>% t
    sw <- val[,sw_ind,drop=FALSE] %>% unlist %>% as.numeric
  }

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

  # Heatmap annotation representing Silhouette width of selected RAVs
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
  )
}
