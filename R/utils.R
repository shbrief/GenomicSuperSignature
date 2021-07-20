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


