context("GenomicSuperSignature functions tests")

data(miniRAVmodel)
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, miniRAVmodel)

test_that("validatedSignatures.R", {
  val_sig <- validatedSignatures(val_all, num.out = 3, scoreCutoff = 0)
  val_sig2 <- validatedSignatures(val_all, num.out = 3, scoreCutoff = 0,
                                  indexOnly = TRUE)

  expect_true(is.matrix(val_sig))
  expect_true(is.vector(val_sig2))
  expect_true(is.numeric(val_sig2))
  expect_equal(val_sig2, c(1076, 2538, 338))
})

expect_error(getModel("C3"), "Prior you entered isn\'t available yet.")
expect_error(meshTable(miniRAVmodel, 1))
expect_equal(dim(meshTable(miniRAVmodel, 1994)), c(31, 2))

test_that("findsignature.R", {
  res <- findSignature(miniRAVmodel, keyword = "Bcell")
  res2 <- findSignature(miniRAVmodel, keyword = "Bcell", k = 5)
  res3 <- findKeywordInRAV(miniRAVmodel, "Bcell", ind = 695)

  expect_true(is.data.frame(res))
  expect_true(is.vector(res2))
  expect_true(is.numeric(res2))
  expect_equal(dim(res), c(2, 2))
  expect_equal(res2, c(16, 17))
  expect_equal(res3, "1|2|3|4|5|6|9")
})

test_that("annotatePCcluster.R", {
  a <- annotateRAV(miniRAVmodel, ind = 695)
  des <- a$Description
  val <- c("IRIS_Bcell-Memory_IgG_IgA", "DMAP_BCELLA3", "IRIS_Bcell-Memory_IgM",
           "IRIS_Bcell-naive", "DMAP_BCELLA4")

  expect_equal(dim(a), c(5, 4))
  expect_equal(colnames(a), c("Description","NES", "pvalue", "qvalues"))
  expect_equal(des, val)
})

test_that("buildAvgLoading.R", {
  data(miniAllZ)
  data(res_hcut)
  res <- buildAvgLoading(miniAllZ, k = 40, cluster = res_hcut$cluster)

  val <- c("cluster", "size", "avgLoading", "k", "n", "studies")

  expect_true(is.list(res))
  expect_equal(names(res), val)
  expect_equal(colnames(res$avgLoading)[40], "Cl40_40 (5/2)")
  expect_equal(names(res$cluster)[1], "GSE13294_eset.PC1")

  val2 <- findStudiesInCluster(miniRAVmodel, 1076)
  expect_equal(length(val2), 10)
  expect_equal(val2[1], "SRP028155")
})
