#' Validation result in data frame
#'
#' @param data An output matrix from \code{\link{validate}} function.
#' @param num.out A number of highly validated PCclusters to output. Default is 5.
#' If \code{cutoff} is provided, this argument is ignored.
#' @param cutoff A numeric value for the minimum correlation. If you don't specify
#' this argument, this function will output the number of highly validated PCclusters
#' defined through \code{num.out}.
#'
#' @return A subset of the input matrix, which meets the given condition.
#' @export
validatedIndex <- function(data, num.out = 5, cutoff = NULL) {
    ind <- which(rownames(data) == "score")
    if (!is.null(cutoff)) {
        validated_ind <- which(data[ind,] > cutoff)
        return(validated_ind)
    } else {
        validated_ind <- order(data[ind,], decreasing = TRUE)[1:num.out]
        return(validated_ind)
    }
}

subsetByCutoff <- function(data, cutoff, subset.by = c("score", "sw")) {
    ind <- which(rownames(data) == subset.by)
    validated_ind <- validatedIndex(data, cutoff = cutoff)
    dat <- data[ind, validated_ind, drop=FALSE]
    return(dat)
}

subsetTop <- function(data, num.out, subset.by = c("score", "sw")) {
    ind <- which(rownames(data) == "score")
    validated_ind <- validatedIndex(data, cutoff = NULL, num.out = num.out)

    ind2 <- which(rownames(data) == subset.by)
    dat <- data[ind2, validated_ind, drop=FALSE]
    return(dat)
}

#' Validation result in heatmap format
#'
#' @param data An output matrix from \code{\link{validate}} function. Subset of
#' this matrix is plotted as a heatmap using \code{\link[ComplexHeatmap]{Heatmap}}
#' @param row_title A character string. Provide the row title.
#' @param num.out A number of highly validated PCclusters to subset to. Default is 5.
#' If \code{cutoff} is provided, this argument is ignored.
#' @param cutoff A numeric value for the minimum correlation. If you don't specify
#' this argument, this function will output the number of highly validated PCclusters
#' defined through \code{num.out}.
#' @param breaks A numeric vector of length 3. Number represents the values assigned
#' to three colors. Default is \code{c(0, 0.5, 1)}.
#' @param colors A character vector of length 3. Each represents the color assigned
#' to three breaks. Default is \code{c("white", "white smoke", "red")}.
#' @param ... any additional argument for \code{ComplexHeatmap::Heatmap}
#'
#' @return A heatmap with the validation result subset through the given conditions.
#'
#' @export
heatmapTable <- function(data, row_title,
                         num.out = 5, cutoff = NULL,
                         breaks = c(0, 0.5, 1),
                         colors = c("white", "white smoke", "red"), ...) {

    # Subset score/silhouette width satisfying the criteria
    if (!is.null(cutoff)) {
        dat <- subsetByCutoff(data, cutoff, subset.by = "score")
        sw <- subsetByCutoff(data, cutoff, subset.by = "sw") %>% as.numeric
        column_title <- paste("Validated with cutoff >", cutoff)
    } else if (is.null(cutoff)) {
        dat <- subsetTop(data, num.out, subset.by = "score")
        sw <- subsetTop(data, num.out, subset.by = "sw") %>% as.numeric
        column_title <- paste("Top", num.out, "validated PCclusters")
    }

    # Define Heatmap Table Size
    if (missing(column_title)) {
        column_title <- character(0)
        hmHeight <- nrow(dat) + 1
    } else {
        column_title <- column_title
        hmHeight <- nrow(dat) + 3
    }

    if (missing(row_title)) {
        row_title <- character(0)
        hmWidth <- ncol(dat) + 1
    } else {
        row_title <- row_title
        hmWidth <- ncol(dat) + 2
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
