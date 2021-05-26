#' Find the RAVs with the keyword-containing enriched pathways
#'
#' This function finds RAVs containing the keyword you provide. If you provide
#' "the number of keyword-containing pathways per RAV" in argument \code{k},
#' it will give you the RAV number.
#'
#' @param RAVmodel PCAGenomicSignatures-object
#' @param keyword A character vector. If you are searching for multiple keywords
#' at the same time, use \code{\link{paste}} with \code{collapse="|"} argument.
#' @param n The number of top ranked (based on abs(NES)) pathways you want to
#' search your keyword
#' @param k The number of keyword-containing pathways you want to get the RAV
#' number. Under default (\code{NULL}), the output will be a data frame with two
#' columns: '# of keyword-containing pathways' and 'Freq'. If you assign the
#' value for this argument, the output will be an integer vector containing the
#' RAV index.
#'
#' @return A data frame or integer vector depending on the parameter \code{k}.
#'
#' @examples
#' data(miniRAVmodel)
#' findSignature(miniRAVmodel, "Bcell")
#' findSignature(miniRAVmodel, "Bcell", k = 5)
#'
#' @export
findSignature <- function(RAVmodel, keyword, n = 5, k = NULL) {
  gsea_all <- gsea(RAVmodel)

  # top n pathways based on abs(NES) value
  topPathways <- sapply(gsea_all, function(x) {
    y <- x[order(abs(x$NES), decreasing = TRUE),][seq_len(n),]
    grep(keyword, y[["Description"]], ignore.case = TRUE)
  })

  nTopPathways <- vapply(topPathways, length, FUN.VALUE = integer(1))
  res <- nTopPathways %>% table %>% as.data.frame
  colnames(res)[1] <- c("# of keyword-containing pathways")

  if (is.null(k)) {
    return(res)
  } else {
    if (!k %in% res[,1]) {
      warning(paste("There is no RAV with", k,
                    "keyword-containing, enriched pathways."))
    } else {
      res <- which(nTopPathways == k) %>% as.numeric
      return(res)
    }
  }
}



#' Find the rank of your keyword in the RAV's GSEA annotation
#'
#' Once you provide RAVmodel, keyword you're searching for, and the RAV number
#' to this function, it will give you the abs(NES)-based rank of your keyword in
#' the enriched pathways of the target RAV. If can be useful to find out how
#' uniquely your keyword-containing pathways are represented.
#'
#' @param RAVmodel PCAGenomicSignatures-object.
#' @param keyword A character vector. If you are searching for multiple keywords
#' at the same time, use \code{\link{paste}} with \code{collapse="|"} argument.
#' @param ind An integer. The RAV number you want to check.
#' @param n An integer. The number of top enriched pathways (based on abs(NES))
#' to search. Under default (\code{NULL}), all the enriched pathways are used.
#' @param includeTotal Under the default condition (\code{TRUE}), the total
#' number of enriched pathways will be also printed out as a part of the output.
#'
#' @return A character containing the rank of keyword-containing pathways
#' (separated by |), followed by the total number of enriched pathways in
#' parenthesis.
#'
#' @examples
#' data(miniRAVmodel)
#' findKeywordInRAV(miniRAVmodel, "Bcell", ind = 695)
#'
#' @export
findKeywordInRAV <- function(RAVmodel, keyword, ind,
                             n = NULL, includeTotal = TRUE) {
 
   .availableRAV(RAVmodel, ind)
  name <- paste0("RAV", ind)
  gsea <- gsea(RAVmodel)[[name]]
  total <- nrow(gsea)

  gsea <- gsea[order(abs(gsea$NES), decreasing = TRUE),,drop = FALSE]
  if (!is.null(n)) {gsea <- gsea[seq_len(n),]}

  keywordRank <- grep(keyword, gsea$Description, ignore.case = TRUE)
  if (length(keywordRank) == 0) {
    warning(paste(name, "doesn't have any pathway with the keyword,", keyword))
  } else {
    if (!includeTotal) {
      res <- paste(keywordRank, collapse = "|")
    } else {
      res <- paste0(paste(keywordRank, collapse = "|"), " (out of ", total, ")")
    }
    return(res)
  }
}
