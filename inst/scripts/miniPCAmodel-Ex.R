## The following is the code used to create this mini PCAmodel from the full
## PCAmodel_PLIERpriors.The full PCAmodel_PLIERpriors was created by the pipeline
## at https://github.com/shbrief/model_building.
\dontrun{
  library(PCAGenomicSignatures)
  library(bcellViper)
  data(bcellViper)
  getModel("PLIERpriors")
  PCAmodel <- readRDS("PCAmodel_PLIERpriors.rds")

  # include top 15 validated RAVs
  val_all <- validate(dset, PCAmodel)
  ind_validated <- validatedSignatures(val_all, num.out = 15, indexOnly = TRUE)

  # include 2 keyword-containing RAVs
  ind_keyword <- findSignature(PCAmodel, "Bcell", k = 5)

  ind_all <- unique(c(ind_validated, ind_keyword))
  miniPCAmodel <- PCAmodel[,ind_all]
  metadata(miniPCAmodel)$size <- metadata(miniPCAmodel)$size[ind_all]

  save(miniPCAmodel, file = "data/miniPCAmodel.RData", compress = "bzip2")
}
