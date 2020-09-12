#' Annotate top PCs from the dataset
#'
#' This function finds the PCcluster scored highest with the top PCs of the dataset,
#' including PCclusters with the negative average silhouette.
#'
#' @param PCnum PC number of your dataset you want to get the annotation results. It
#' should be an integer value between 1 and 8.
#' @param val_all The output from \code{\link{validate}}
#' @param PCAmodel The PCAmodel used to generate the argument, \code{val}
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
#'
#' @return A data frame of a list based on the \code{simplify} argument. Check the
#' output detail above.
#'
#' @export
annotatePC <- function(PCnum, val_all, PCAmodel, n = 5,
                       scoreCutoff = 0.5, nesCutoff = NULL,
                       simplify = TRUE) {
  res <- vector(mode = "list", length = length(PCnum))

  for (i in seq_along(PCnum)) {
    annotPC <- validatedSignatures(val_all, num.out = 1, scoreCutoff = scoreCutoff, whichPC = PCnum[i])
    cl_name <- paste0("PCcluster", annotPC[,"cl_num"])
    annotatedCluster <- gsea(PCAmodel)[[cl_name]]

    if (is.null(annotatedCluster)) {
      emptyTable <- gsea(PCAmodel)[[1]][0,] # assign the empty gsea table
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

      topAnnotation <- topAnnotation[1:n,]
      rownames(topAnnotation) <- NULL
      res[[i]] <- topAnnotation
      names(res)[i] <- paste0("PC", PCnum[i], "-", cl_name)
    }
  }

  if (isTRUE(simplify)) {
    simple_res <- sapply(res, function(x) x$Description) %>% as.data.frame
    return(simple_res)
  } else {
    return(res)
  }

}


annotatePC_internal <- function(PCnum, dat, n = 5, simplify = TRUE, gseaDB = NULL) {
  a <- readRDS(file.path(dat_dir, gseaDB))
  res <- vector(mode = "list", length = length(PCnum))

  for (i in seq_along(PCnum)) {
    annotPC <- validatedSignatures(dat, num.out = 1, whichPC = PCnum[i])
    cl_name <- paste0("PCcluster_", annotPC[,"cl_num"])
    annotatedCluster <- a[[cl_name]]
    topAnnotation <- annotatedCluster[order(abs(annotatedCluster$NES), decreasing = TRUE),][1:n,]
    rownames(topAnnotation) <- NULL
    res[[i]] <- topAnnotation
    names(res)[i] <- paste0("PC", PCnum[i], "-", cl_name)
  }

  if (isTRUE(simplify)) {
    simple_res <- sapply(res, function(x) x$Description)
    return(simple_res)
  } else {
    return(res)
  }
}
