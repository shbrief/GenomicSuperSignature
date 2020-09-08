#' Two-dimentional PCA plot with the PC annotation
#'
#' @param dataset An expression matrix with genes (rows) x samples (columns)
#' @param PCAmodel PCAGenomicSignatures-class object
#' @param PCs A numeric vector length of 2. It should be between 1 and 8.
#' @param val_all The output from \code{\link{validate}}
#' @param color_by A named vector with the feature you want to color by. Name should
#' be match with the sample names of the dataset.
#' @param color_lab A name for color legend. If this argument is not provided, the
#' color legend will be labeld as "Color By" by default.
#' @param trimed_pathway_len Positive inter values, which is the display width of
#' pathway names. Default is 45.
#'
#' @importFrom ggpubr ggtexttable ttheme table_cell_font ggarrange tab_nrow
#' @importFrom ggplot2 ggplot geom_point
#'
#' @export
plotAnnotatedPCA <- function(dataset, PCAmodel, PCs, val_all = NULL, color_by = NULL, color_lab = NULL,
                            trimed_pathway_len = 45) {
  if (is.null(val_all)) {val_all <- validate(dataset, PCAmodel)}
  PCAres <- prcomp(dataset)$rotation %>% as.data.frame

  # two PCs to plot
  ind1 <- which(colnames(PCAres) == paste0("PC", PCs[1]))
  ind2 <- which(colnames(PCAres) == paste0("PC", PCs[2]))

  if (is.null(color_lab)) {color_lab = "Color By"}
  if (!is.null(color_by)) {
    colorFeature <- stack(color_by) %>% tibble::column_to_rownames(var = "ind")
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
  annotatedPC <- annotatePC(c(ind1, ind2), val_all, PCAmodel)
  a <- which(nchar(annotatedPC[,1]) > trimed_pathway_len)
  annotatedPC[a,1] <- paste0(strtrim(annotatedPC[a,1], trimed_pathway_len), "...")
  b <- which(nchar(annotatedPC[,2]) > trimed_pathway_len)
  annotatedPC[b,2] <- paste0(strtrim(annotatedPC[b,2], trimed_pathway_len), "...")

  # PC annotation table
  myTable <- ggpubr::ggtexttable(annotatedPC,
                                 rows = NULL, theme = ggpubr::ttheme("mOrange")) %>%
    ggpubr::table_cell_font(row = 2:ggpubr::tab_nrow(.), column = 1:2, size = 9)

  res <- ggpubr::ggarrange(myPlot, myTable, ncol = 1, nrow = 2, heights = c(1, 0.5))
  print(res)
}

