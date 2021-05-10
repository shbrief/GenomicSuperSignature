#' Common genes from datasets
#'
#' @param dat A character vector containing the name of datasets loaded in memory.
#' Currently, the loaded dataset format should be in ExpressionSet class. You
#' can also provide a list of genes you want to merge into a common gene list.
#' @param loaded The default condition (\code{TRUE}) assumes that the dataset (e.g.
#' expression matrix with genes as row names) is loaded. Set this argument as
#' \code{FALSE} if you provide a list of genes in the argument \code{dat}.
#'
#' @return A character vector with the common genes of the datasets.
#' @note This function is for CRC datasets.
#'
commonGene <- function(dat, loaded = TRUE) {
    if (loaded) {
        lapply(dat, function(x) {
            dat <- get(x)
            rownames(Biobase::exprs(dat))
        })
        Reduce(intersect, dat)
    } else {
        Reduce(intersect, dat)
    }
}

