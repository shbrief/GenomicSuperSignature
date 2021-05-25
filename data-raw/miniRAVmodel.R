## The following is the code used to create this miniRAVmodel from the full
## RAVmodel_PLIERpriors.The full RAVmodel_PLIERpriors was created by the
## pipeline at https://github.com/shbrief/model_building.
\dontrun{
  library(GenomicSuperSignature)
  library(bcellViper)
  data(bcellViper)
  getModel("PLIERpriors")
  RAVmodel <- readRDS("RAVmodel_PLIERpriors.rds")

  ## RAVs for miniRAVmodel
  # top 15 validated RAVs
  val_all <- validate(dset, RAVmodel)
  ind_validated <- validatedSignatures(val_all, num.out=15, indexOnly=TRUE)
  # two keyword-containing RAVs
  ind_keyword <- findSignature(RAVmodel, "Bcell", k=5)
  ind_all <- unique(c(ind_validated, ind_keyword))

  ## Subset RAVmodel
  # subset RAVindex
  miniRAVmodel <- RAVmodel[,ind_all]
  # subset metadata$cluster
  cl <-  which(metadata(miniRAVmodel)$cluster %in% ind_all)
  metadata(miniRAVmodel)$cluster <- metadata(miniRAVmodel)$cluster[cl]
  # subset metadata$size
  metadata(miniRAVmodel)$size <- metadata(miniRAVmodel)$size[ind_all]

  # check the best compression using `tools::checkRdaFiles`
  save(miniRAVmodel, file = "data/miniRAVmodel.RData", compress = "bzip2")
}
