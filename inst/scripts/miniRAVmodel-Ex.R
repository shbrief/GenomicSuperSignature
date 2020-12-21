## The following is the code used to create this miniRAVmodel from the full
## RAVmodel_PLIERpriors.The full RAVmodel_PLIERpriors was created by the pipeline
## at https://github.com/shbrief/model_building.
\dontrun{
  library(GenomicSuperSignature)
  library(bcellViper)
  data(bcellViper)
  RAVmodel <- getModel("PLIERpriors")

  # include top 15 validated RAVs
  val_all <- validate(dset, RAVmodel)
  ind_validated <- validatedSignatures(val_all, num.out = 15, indexOnly = TRUE)

  # include 2 keyword-containing RAVs
  ind_keyword <- findSignature(RAVmodel, "Bcell", k = 5)

  ind_all <- unique(c(ind_validated, ind_keyword))
  miniRAVmodel <- RAVmodel[,ind_all]
  metadata(miniRAVmodel)$size <- metadata(miniRAVmodel)$size[ind_all]

  save(miniRAVmodel, file = "data/miniRAVmodel.RData", compress = "bzip2")
}
