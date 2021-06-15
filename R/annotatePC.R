#' Annotate top PCs from the dataset
#'
#' This function finds the RAV with the highest validation score (including
#' RAVs with negative silhouette width) for specified PC of the dataset and
#' returns the top enriched pathways.
#'
#' @param PCnum A numeric vector. PC number of your dataset that you want to get
#' the annotation results. The vector can contain any integer number among
#' \code{1:8}.
#' @param val_all The output from \code{\link{validate}}
#' @param RAVmodel The RAVmodel used to generate the input for the argument,
#' \code{val_all}.
#' @param n An integer. Default is 5. The number of the top enriched pathways
#' to print out. If there are fewer than n pathways passed the cutoff, it will
#' print out \code{NA}.
#' @param scoreCutoff A numeric value for the minimum correlation.
#' Default is 0.5.
#' @param nesCutoff A numeric value for the minimum NES. Default is \code{NULL}
#' and the suggested value is 3.
#' @param simplify A logical. Under default (\code{TRUE}), the output will be a
#' data frame with the number of column same as the length of \code{PCnum}
#' argument, and the number of row same as the \code{n} argument. If it is set
#' to \code{FALSE}, the output will be a list with the length of \code{PCnum}
#' argument, where each element is a data frame containing detailed GSEA output
#' of enriched pathways.
#' @param abs Default is \code{FALSE}. If it's set to \code{TRUE}, the enriched
#' pathways will be listed based on \code{abs(NES)}.
#' @param trimed_pathway_len Positive inter values, which is the display width
#' of pathway names. Default is 45.
#'
#' @return A data frame of a list based on the \code{simplify} argument. Check
#' the output detail above.
#'
#' @examples
#' data(miniRAVmodel)
#' library(bcellViper)
#' data(bcellViper)
#' val_all <- validate(dset, miniRAVmodel)
#' annotatePC(2, val_all, miniRAVmodel)
#'
#' @export
annotatePC <- function(PCnum, val_all, RAVmodel, n = 5,
                       scoreCutoff = 0.5, nesCutoff = NULL,
                       simplify = TRUE, abs = FALSE,
                       trimed_pathway_len = 45) {

  ## Input validation
  if(any(!PCnum %in% seq(8))) {stop("PCnum should be an integer among c(1:8).")}
  stopifnot(length(simplify) == 1L, !is.na(simplify), is.logical(simplify))
  stopifnot(length(abs) == 1L, !is.na(abs), is.logical(abs))

  res <- vector(mode = "list", length = length(PCnum))

  for (i in seq_along(PCnum)) {
    annotPC <- validatedSignatures(val_all, num.out = 1,
                                   scoreCutoff = scoreCutoff,
                                   whichPC = PCnum[i])
    cl_name <- paste0("RAV", annotPC[,"cl_num"])
    annotatedCluster <- gsea(RAVmodel)[[cl_name]]

    if (is.null(annotatedCluster)) {
      emptyTable <- gsea(RAVmodel)[[1]][0,] # assign the empty GSEA table
      emptyTable[seq(n),1] <- "No significant pathways"
      res[[i]] <- emptyTable
      names(res)[i] <- paste0("PC", PCnum[i], "-noAnnot")
    } else {
      # absolute value of NES
      FUN <- function(x) {if (abs) {abs(x)} else {I(x)}}

      # apply NES cutoff
      if (!is.null(nesCutoff)) {
        annotatedCluster <- annotatedCluster[
          FUN(annotatedCluster$NES) >= nesCutoff,,drop = FALSE]
      }

      # subset GSEA results
      topAnnotation <- annotatedCluster[
        order(FUN(annotatedCluster$NES), decreasing = TRUE),,drop = FALSE]
      topAnnotation <- topAnnotation[seq_len(n),]
      rownames(topAnnotation) <- NULL

      res[[i]] <- topAnnotation
      names(res)[i] <- paste0("PC", PCnum[i], "-", cl_name)
    }
  }

  # Trim the long pathway names
  for (i in seq_along(res)) {
    ind <- which(nchar(res[[i]]$Description) > trimed_pathway_len)
    res[[i]]$Description[ind] <- paste0(strtrim(res[[i]]$Description[ind],
                                                trimed_pathway_len), "...")
  }

  # Output formats:
  # A list of data frames with details vs. data frame with description only
  if (simplify) {
    simple_res <- lapply(res, function(x) x$Description) %>% as.data.frame
    return(simple_res)
  } else {
    return(res)
  }
}
