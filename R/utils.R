### Extract expression matrix from different classes of input datasets
.extractExprsMatrix <- function(dataset) {
  if (is(dataset, "ExpressionSet")) {
    dat <- Biobase::exprs(dataset)
  } else if (is(dataset,"SummarizedExperiment")) {
    dat <- SummarizedExperiment::assay(dataset)
  } else if (is.matrix(dataset)) {
    dat <- dataset
  } else {
    stop("'dataset' should be one of the following objects: ExpressionSet,
         SummarizedExperiment, and matrix.")
  }
  return(dat)
}


### Check ind validity
.availableRAV <- function(RAVmodel, ind) {
  availableRAV <- gsub("RAV", "", colData(RAVmodel)$RAV) %>% as.numeric

  ## Check whether the ind exists in the model
  x <- vector(length = length(ind))
  for (i in seq_along(ind)) {
    if (!ind[i] %in% availableRAV) {
      x[i] <- TRUE   # assign TRUE if index doesn't exist
    }
  }

  ## Print error message if any of the ind doesn't exist.
  if (any(x)) {
    y <- paste(paste0("RAV",ind[x]), collapse=", ")  # combine non-existing ind
    msg <- paste0("Selected ind (", y, ") doesn't exist.")
    stop(msg)
  }
}


## Restructure RAVmodel metadata slot
.RAVmodelVersion <- function(RAVmodel) {
  if (version(RAVmodel) == ">= 0.0.7") {
    cluster <- S4Vectors::metadata(RAVmodel)$cluster
  } else {
    cluster <- colData(RAVmodel)$cluster
  }
}


## Extract variance explained by PCs in a given cluster
.varByPCsInCluster <- function(RAVmodel, ind) {
  # components in clusters
  cl_membership <- metadata(RAVmodel)$cluster
  components <- names(which(cl_membership == ind))

  # PCA summary
  pcaSummary <- trainingData(RAVmodel)$PCAsummary
  Projs <- lapply(components, function(x) {
    unlist(strsplit(x, "\\.PC"))[1] %>% as.character
  }) %>% unlist
  data <- pcaSummary[Projs]

  # Extract variance explained
  input_summary <- as.data.frame(matrix(ncol = 3, nrow = length(data)))
  colnames(input_summary) <- c("studyName", "PC", "Variance explained (%)")

  for (i in seq_along(data)) {
    studyname <- Projs[i]
    j <- unlist(strsplit(components[i], "\\.PC"))[2] %>% as.numeric
    var <- data[[i]]["Variance",j]

    input_summary[i, 1] <- studyname
    input_summary[i, 2] <- j
    input_summary[i, 3] <- round(var*100, digits = 2)
  }

  return(input_summary)
}


## Message for low-quality RAVs
.lowQualityRAVs <- function(RAVmodel, ind, filterMessage = TRUE) {

  if (isTRUE(filterMessage)) {
    ## Load filterList
    local_data_store <- new.env(parent = emptyenv())
    data("filterList", envir = local_data_store, package = "GenomicSuperSignature")
    filterList <- local_data_store[["filterList"]]

    ## Select RAVmodel
    filterListNames <- c("Cluster_Size_filter", "GSEA_C2_filter",
                         "GSEA_PLIERpriors_filter", "Redundancy_filter")
    c2 <- "MSigDB C2 version 7.1"
    plier_priors <- "Three priors from PLIER (bloodCellMarkersIRISDMAP, svmMarkers, and canonicalPathways)"

    if (nrow(trainingData(RAVmodel)) == 536 & geneSets(RAVmodel) == c2) {
      filterList <- filterList[filterListNames[c(1,2,4)]]
    } else if ((nrow(trainingData(RAVmodel)) == 536 &
                geneSets(RAVmodel) == plier_priors)) {
      filterList <- filterList[filterListNames[c(1,3,4)]]
    }

    ## Check whether index belong to the filter list
    for (i in ind) {
      res <- vapply(filterList, function(x) {i %in% x}, logical(1))
      if (any(res)) {
        filtered <- paste(names(res)[which(res == TRUE)], collapse = ", ") %>%
          gsub("_filter", "", .)
        msg <- paste(paste0("RAV", i), "can be filtered based on", filtered)
        message(msg)
      }
    }

    ## More information on GenomicSuperSignaturePaper GitHub page
    # if (any(res)) {message("Information on filtering : bit.ly/rav_filtering")}
  }
}


## Study metadata for different RAVmodels
.getStudyMeta <- function(RAVmodel) {

  td <- rownames(trainingData(RAVmodel)) # training data used for RAVmodel
  dir <- system.file("extdata", package = "GenomicSuperSignature")

  if ("DRP000987" %in% td) {
    ## 536 datasets from refine.bio
    studyMeta <- utils::read.table(file.path(dir, "studyMeta.tsv.gz"))
  } else if ("GSE13294" %in% td) {
    ## 8 CRC and 10 OV from curated data packages
    studyMeta <- utils::read.table(file.path(dir, "studyMeta_CRCOV.tsv"))
  }

  return(studyMeta)
}
