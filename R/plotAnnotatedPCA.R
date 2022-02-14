#' Two-dimentional PCA plot with the PC annotation
#'
#' @import ggpubr
#' @import ggplot2
#'
#' @param dataset A gene expression profile to be validated. Different classes
#' of objects can be used including ExpressionSet, SummarizedExperiment,
#' RangedSummarizedExperiment, or matrix. Rownames (genes) should be in symbol
#' format. If it is a matrix, genes should be in rows and samples in columns.
#' @param RAVmodel PCAGenomicSignatures-class object
#' @param PCnum A numeric vector length of 2. It should be between 1 and 8.
#' @param val_all The output from \code{\link{validate}}
#' @param scoreCutoff A numeric value for the minimum correlation.
#' Default 0.5.
#' @param nesCutoff A numeric value for the minimum NES. Default is \code{NULL}
#' and the suggested value is 3.
#' @param color_by A named vector with the feature you want to color by. Name
#' should be match with the sample names of the dataset.
#' @param color_lab A name for color legend. If this argument is not provided,
#' the color legend will be labeled as "Color By" by default.
#' @param trimed_pathway_len Positive inter values, which is the display width
#' of pathway names. Default is 45.
#'
#' @return Scatter plot and the table with annotation. If enriched pathway
#' didn't pass the \code{scoreCutoff} the table will be labeled as "No
#' significant pathways". If any enriched pathway didn't pass the
#' \code{nesCutoff}, it will labeled as NA.
#'
#' @examples
#' data(miniRAVmodel)
#' library(bcellViper)
#' data(bcellViper)
#' \dontrun{
#' plotAnnotatedPCA(exprs(dset), miniRAVmodel, PCnum = c(1,2))
#' }
#'
#' @export
plotAnnotatedPCA <- function(dataset, RAVmodel, PCnum, val_all = NULL,
                             scoreCutoff = 0.5, nesCutoff = NULL,
                             color_by = NULL, color_lab = NULL,
                             trimed_pathway_len = 45) {

    # Extract expression matrix from different classes
    dat <- .extractExprsMatrix(dataset)

    if (is.null(val_all)) {val_all <- validate(dat, RAVmodel)}
    PCAres <- stats::prcomp(dat)$rotation %>% as.data.frame

    # two PCs to plot
    ind1 <- which(colnames(PCAres) == paste0("PC", PCnum[1]))
    ind2 <- which(colnames(PCAres) == paste0("PC", PCnum[2]))

    if (is.null(color_lab)) {color_lab <- "Color By"}
    if (!is.null(color_by)) {
        colorFeature <- utils::stack(color_by) %>%
            tibble::column_to_rownames(var = "ind")
        PCAres <- cbind(PCAres, colorFeature)
        myPlot <- ggplot(PCAres, mapping = aes_string(x = names(PCAres)[ind1],
                                                      y = names(PCAres)[ind2],
                                                      color = "values")) +
            labs(color = color_lab) +
            geom_point()
    } else {
        myPlot <- ggplot(PCAres, mapping = aes_string(x = names(PCAres)[ind1],
                                                      y = names(PCAres)[ind2])) +
            geom_point()
    }

    # Trim the long pathway names
    annotatedPC <- annotatePC(c(ind1, ind2), val_all, RAVmodel,
                              scoreCutoff = scoreCutoff, nesCutoff = nesCutoff,
                              trimed_pathway_len = trimed_pathway_len)

    # PC annotation table - flextable version
    flextable::set_flextable_defaults(na_str = "NA")
    myTable <- flextable::flextable(annotatedPC) %>% flextable::width(width = 2.5)
    myTable <- grid::rasterGrob(flextable::as_raster(myTable))

    # # PC annotation table - ggpubr version
    # myTable <- ggpubr::ggtexttable(annotatedPC,
    #                                rows = NULL,
    #                                theme = ggpubr::ttheme("mOrange"))
    # myTable <- ggpubr::table_cell_font(myTable,
    #                                    row = 2:ggpubr::tab_nrow(myTable),
    #                                    column = 1:2, size = 9)

    res <- ggpubr::ggarrange(myPlot, myTable,
                             ncol = 1, nrow = 2, heights = c(1, 0.5))
    print(res)
}

