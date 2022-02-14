context("Test annotateRAV function")

data(miniRAVmodel)

test_that("Confirm the output structure", {
  res <- annotateRAV(miniRAVmodel, ind = 695)
  des <- res$Description
  ## RAVmodel version before 1.1.1
  # val <- c("IRIS_Bcell-Memory_IgG_IgA", "DMAP_BCELLA3",
  #          "IRIS_Bcell-Memory_IgM", "IRIS_Bcell-naive", "DMAP_BCELLA4")
  ## RAVmodel version above 1.1.1
  val <- c("IRIS_Bcell-Memory_IgG_IgA", "DMAP_BCELLA3",
           "IRIS_Bcell-naive", "IRIS_Bcell-Memory_IgM", "DMAP_BCELLA2")

  expect_equal(dim(res), c(5, 4))
  expect_equal(colnames(res), c("Description","NES", "pvalue", "qvalues"))
  expect_equal(des, val)
})

test_that("Test abs argument", {
  res2 <- annotateRAV(miniRAVmodel, ind = 1076)
  res3 <- annotateRAV(miniRAVmodel, ind = 1076, abs = TRUE)
  expect_equal(res2[1,1], "REACTOME_METABOLISM_OF_PROTEINS")
  expect_equal(res3[1,1], "REACTOME_CELL_CYCLE")
})
