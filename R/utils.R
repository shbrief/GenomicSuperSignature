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
    unlist(strsplit(x, "\\."))[1] %>% as.character
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
