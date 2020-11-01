#' Select enriched pathways for PCAmodel
#'
#' From the GSEA output of each RAV, subset the enriched pathways with the
#' minimum q-value. Also, this function keeps only \code{Description, NES, qvalues}.
#'
#' @param PCAmodel PCAGenomicSignatures object.
#' @param gsea.dir "~/data2/PCAGenomicSignatureLibrary/refinebioRseq/PCAmodel_536/gsea"
#' @return A list with the length of RAVs in the provided PCAmodel. Each
#' element contains the subset data.frame of GSEA output.
#'
#' @note This function is for model construction, not for end-users.
#'
searchPathways <- function(PCAmodel, gsea.dir) {

  ## If you want to select only a subset of RAVs with the specific cluster size
  # ind <- which(metadata(PCAmodel)$size > 3)
  # gsea_all <- vector(mode = "list", length = length(ind))
  # names(gsea_all) <- colnames(PCAmodel)[ind]

  gsea_all <- vector(mode = "list", length = ncol(PCAmodel))
  names(gsea_all) <- colnames(PCAmodel)
  gsea.dir <- gsea.dir

  for (i in seq_len(ncol(PCAmodel))) {
    pathToRes <- file.path(gsea.dir, paste0("gsea_", i, ".rds"))
    res <- readRDS(pathToRes)

    # If there is no enriched pathways
    if (nrow(res) == 0) {
      resName <- paste0("RAV", i)
      gsea_all[[resName]] <- NA
      print(paste("RAV", i, "has no enriched pathways."))
      next
    }

    res <- res[which(res$qvalues == min(res$qvalues)), c("Description", "NES", "pvalue", "qvalues"), drop = FALSE]
    resName <- paste0("RAV", i)
    gsea_all[[resName]] <- res
    print(paste("RAV", i, "is added."))
  }

  return(gsea_all)
}
