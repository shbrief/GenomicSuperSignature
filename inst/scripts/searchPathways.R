#' Select enriched pathways for RAVmodel
#'
#' From the GSEA output of each RAV, subset the enriched pathways with the
#' minimum q-value. Also, this function keeps only \code{Description, NES, qvalues}.
#'
#' @param RAVmodel PCAGenomicSignatures object
#' @param gsea.dir "~/data2/PCAGenomicSignatureLibrary/refinebioRseq/RAVmodel_536/gsea"
#' @param min.qval Under the default (\code{TRUE}), only the enriched pathways
#' with the minimum qvalue will be saved. If it is set \code{FALSE}, all the
#' enriched pathways pass the \code{pvalueCutoff} of the function
#' \code{\link[clusterProfiler]{GSEA}} will be extracted.
#'
#' @return A list with the length of RAVs in the provided RAVmodel. Each
#' element contains the subset data.frame of GSEA output.
#'
#' @note This function is for model construction, not for end-users.
#'
searchPathways <- function(RAVmodel, gsea.dir, min.qval = TRUE) {

  ## If you want to select only a subset of RAVs with the specific cluster size
  # ind <- which(metadata(RAVmodel)$size > 3)
  # gsea_all <- vector(mode = "list", length = length(ind))
  # names(gsea_all) <- colnames(RAVmodel)[ind]

  gsea_all <- vector(mode = "list", length = ncol(RAVmodel))
  names(gsea_all) <- paste0("RAV", seq_len(ncol(RAVmodel)))
  gsea.dir <- gsea.dir

  for (i in seq_len(ncol(RAVmodel))) {
    pathToRes <- file.path(gsea.dir, paste0("gsea_", i, ".rds"))
    res_ls <- readRDS(pathToRes)
    res <- res_ls@result

    ## If there is no enriched pathways
    if (nrow(res) == 0) {
      resName <- paste0("RAV", i)
      gsea_all[[resName]] <- res[, c("Description", "NES", "pvalue", "qvalue"), drop = FALSE]
      print(paste("RAV", i, "has no enriched pathways."))
      next
    }

    ## Select all enriched pathways or only the ones with the minimum qvalue
    if (isTRUE(min.qval)) {
        
        res <- res[which(res$qvalue == min(res$qvalue)),
                   c("Description", "NES", "pvalue", "qvalue"), drop = FALSE]
    } else {
        res <- res[, c("Description", "NES", "pvalue", "qvalue"), drop = FALSE]
    }

    resName <- paste0("RAV", i)
    gsea_all[[resName]] <- res
    print(paste("RAV", i, "is added."))
  }

  return(gsea_all)
}
