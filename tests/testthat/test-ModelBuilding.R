context("Test functions for model buildling")


test_that("getModel function", {
  expect_error(getModel("C3"), "Prior you entered isn\'t available yet.")
})

test_that("buildAvgLoading function", {
  data(miniAllZ)
  data(res_hcut)
  res <- buildAvgLoading(miniAllZ, k = 40, cluster = res_hcut$cluster)
  val <- c("cluster", "size", "avgLoading", "k", "n", "studies")

  expect_true(is.list(res))
  expect_equal(names(res), val)
  expect_equal(colnames(res$avgLoading)[40], "Cl40_40 (5/2)")
  expect_equal(names(res$cluster)[1], "GSE13294_eset.PC1")
})


test_that("findStudiesInCluster function in buildAvgLoading.R", {
  data(miniRAVmodel)
  val2 <- findStudiesInCluster(miniRAVmodel, 1076)
  
  expect_equal(length(val2), 10)
  expect_equal(val2[1], "SRP028155")
})

