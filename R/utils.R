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
  if (!ind %in% availableRAV) {stop("Selected ind doesn't exist.")}
}


