#' Common genes from datasets
#'
#' @import magrittr
#'
#' @param dat A character vector with name of the datasets loaded in memory. You
#' can also provide a list of genes you want to merge into a common gene list.
#' @param loaded Under the default condition (\code{TRUE}), the dataset (e.g. expression
#' matrix with genes as row names) is loaded. \code{FALSE} if you provide a list
#' of gene in the argument \code{dat}.
#'
#' @return A character vector with the common genes of the datasets.
#' @export
commonGene <- function(dat, loaded = TRUE) {
    if (isTRUE(loaded)) {
        lapply(dat, function(x) {
            dat <- get(x)
            rownames(exprs(dat))
        }) %>% Reduce(intersect, .)
    } else if (isFALSE(loaded)) {
        Reduce(intersect, dat)
    }
}

