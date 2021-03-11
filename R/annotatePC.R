#' Annotate top PCs from the dataset
#'
#' This function finds the RAV scored highest with the top PCs of the dataset,
#' including RAVs with the negative average silhouette.
#'
#' @param PCnum A numeric vector. PC number of your dataset that you want to get
#' the annotation results. The vector can contain any integer number among \code{1:8}.
#' @param val_all The output from \code{\link{validate}}
#' @param RAVmodel The RAVmodel used to generate the input for the argument, \code{val_all}.
#' @param n An integer. Default is 5. The number of the top enriched pathaways to
#' print out. If there are fewer than n pathways passed the cutoff, it will print out \code{NA}.
#' @param scoreCutoff A numeric value for the minimum correlation. Default is 0.5.
#' @param nesCutoff A numeric value for the minimum NES. Default is \code{NULL} and
#' the suggested value is 3.
#' @param simplify A logical. Under default (\code{TRUE}), the output will be a
#' data frame with the number of column same as the length of \code{PCnum} argument,
#' and the number of row same as the \code{n} argument. If it is set to \code{FALSE},
#' the output will be a list with the length of \code{PCnum} argument, where each
#' element is a data frame containing detailed GSEA output of enriched pathways.
#' @param abs Default is \code{FALSE}. If it's set to \code{TRUE}, the enriched
#' pathways will be listed based on \code{abs(NES)}.
#' @param trimed_pathway_len Positive inter values, which is the display width of
#' pathway names. Default is 45.
#'
#' @return A data frame of a list based on the \code{simplify} argument. Check the
#' output detail above.
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

  ## Check PCnum is a valid input
  if (any(!PCnum %in% 1:8)) {stop("PCnum should be an integer among c(1:8).")}

  res <- vector(mode = "list", length = length(PCnum))

  for (i in seq_along(PCnum)) {
    annotPC <- validatedSignatures(val_all, num.out = 1, scoreCutoff = scoreCutoff, whichPC = PCnum[i])
    cl_name <- paste0("RAV", annotPC[,"cl_num"])
    annotatedCluster <- gsea(RAVmodel)[[cl_name]]

    if (is.null(annotatedCluster)) {
      emptyTable <- gsea(RAVmodel)[[1]][0,] # assign the empty gsea table
      emptyTable[1:5,1] <- "No significant pathways"
      res[[i]] <- emptyTable
      names(res)[i] <- paste0("PC", PCnum[i], "-noAnnot")
    } else {
      # topAnnotation <- annotatedCluster[order(abs(annotatedCluster$NES), decreasing = TRUE),,drop = FALSE][1:n,]
      if (!is.null(nesCutoff)) {
        topAnnotation <- annotatedCluster[annotatedCluster$NES >= nesCutoff,,drop = FALSE]
        topAnnotation <- topAnnotation[order(topAnnotation$NES, decreasing = TRUE),,drop = FALSE]
      } else {
        topAnnotation <- annotatedCluster[order(annotatedCluster$NES, decreasing = TRUE),,drop = FALSE]
      }

      if (isTRUE(abs)) {
        topAnnotation <- topAnnotation[order(abs(annotatedCluster$NES), decreasing = TRUE),,drop = FALSE]
      }

      topAnnotation <- topAnnotation[1:n,]
      rownames(topAnnotation) <- NULL
      res[[i]] <- topAnnotation
      names(res)[i] <- paste0("PC", PCnum[i], "-", cl_name)
    }
  }

  # Trim the long pathway names
  for (i in seq_along(res)) {
    ind <- which(nchar(res[[i]]$Description) > trimed_pathway_len)
    res[[i]]$Description[ind] <- paste0(strtrim(res[[i]]$Description[ind], trimed_pathway_len), "...")
  }

  # Output format: list of data frame with detail vs. data frame only with description
  if (isTRUE(simplify)) {
    simple_res <- sapply(res, function(x) x$Description) %>% as.data.frame
    return(simple_res)
  } else {
    return(res)
  }
}


# annotatePC_internal <- function(PCnum, dat, n = 5, simplify = TRUE, gseaDB = NULL) {
#   a <- readRDS(file.path(dat_dir, gseaDB))
#   res <- vector(mode = "list", length = length(PCnum))
#
#   for (i in seq_along(PCnum)) {
#     annotPC <- validatedSignatures(dat, num.out = 1, whichPC = PCnum[i])
#     cl_name <- paste0("RAV_", annotPC[,"cl_num"])
#     annotatedCluster <- a[[cl_name]]
#     topAnnotation <- annotatedCluster[order(abs(annotatedCluster$NES), decreasing = TRUE),][1:n,]
#     rownames(topAnnotation) <- NULL
#     res[[i]] <- topAnnotation
#     names(res)[i] <- paste0("PC", PCnum[i], "-", cl_name)
#   }
#
#   if (isTRUE(simplify)) {
#     simple_res <- sapply(res, function(x) x$Description)
#     return(simple_res)
#   } else {
#     return(res)
#   }
# }
