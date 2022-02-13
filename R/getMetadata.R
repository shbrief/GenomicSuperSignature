#' Extract information on a specific RAV
#'
#' @param RAVmodel A PCAGenomicSignatures object
#' @param ind An index of RAV
#' @return A list with four elements: clusterSize, silhouetteWidth,
#' enrichedPathways (the number of enriched pathways), and members.
#' The 'members' is the summary table of PCs in RAV, containing three columns:
#' studyName, PC, and Variance explained (%).
#'
#' @examples
#' data(miniRAVmodel)
#' getRAVInfo(miniRAVmodel, ind = 438)
#'
#' @export
getRAVInfo <- function(RAVmodel, ind) {

    rav <- paste0("RAV", ind)
    meta <- metadata(RAVmodel)
    colDat <- colData(RAVmodel)
    sw <- colDat$silhouetteWidth

    # PCs in the cluster
    PCmembership <- findStudiesInCluster(RAVmodel, ind)

    # metadata with study
    res <- vector(mode = "list")
    res$clusterSize <- as.numeric(meta$size[rav])
    res$silhouetteWidth <- round(as.numeric(sw[as.character(ind)]), 2)
    res$enrichedPathways <- nrow(colDat$gsea[[rav]])
    res$members <- PCmembership

    return(res)
}

#' Extract information on a specific training dataset
#'
#' @param RAVmodel A PCAGenomicSignatures object
#' @param study A character for SRA study accession.
#' @return A list with three elements: studyTitle, studySize (the number of
#' samples from this study used in the RAVmodel building), and RAVs. 'RAVs' is
#' a data frame with three columns - PC (1 to 20), RAV (RAV that the given PC
#' belongs to), and Variance explained (%). In the example, we used
#' miniRAVmodel, which doesn't have all the PCA summary information, so the
#' example will return only the two PCs of the study instead of all twenty.
#'
#' @examples
#' data(miniRAVmodel)
#' getStudyInfo(miniRAVmodel, "SRP028155")
#'
#' @export
getStudyInfo <- function(RAVmodel, study) {
    meta <- metadata(RAVmodel)
    trainingDat <- trainingData(RAVmodel)

    # Index of the requested study
    x <- gsub("\\.PC.*$", "", names(meta$cluster))
    ind <- which(x == study)

    # Variance explained by each PCs
    pcaRes <- trainingDat$PCAsummary[[study]]
    varExplained <- round(pcaRes["Variance",]*100, 2)

    # Clusters where 20 PCs are belong to
    metaCl <- meta$cluster[ind]
    PCnum <- as.numeric(gsub(".*\\.PC", "", names(metaCl)))

    PCdistribution <- data.frame(PC = PCnum,
                                 RAV = as.numeric(metaCl),
                                 Var = varExplained[PCnum],
                                 row.names = NULL)
    colnames(PCdistribution)[3] <- "Variance explained (%)"

    # Study info
    studyMeta <- .getStudyMeta(RAVmodel)
    studyMeta_ind <- which(studyMeta$studyName == study)
    studyTitle <- studyMeta$title[studyMeta_ind]  # title
    studySize <- studyMeta$downloaded[studyMeta_ind]  # size

    # Output
    res <- vector(mode = "list")
    res$studyTitle <- studyTitle
    res$studySize <- studySize
    res$RAVs <- PCdistribution
    return(res)
}
