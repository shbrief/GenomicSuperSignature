context("Test validate function handles differnet types of inputs")

data(miniRAVmodel)
data(miniTCGA)
library(bcellViper)
data(bcellViper)
microTCGA <- readRDS("microTCGA.rds")


test_that("validate works with a single dataset", {
  val_all <- validate(dset, miniRAVmodel)
  val_se <- validate(miniTCGA[[1]], miniRAVmodel)

  expect_true(is.data.frame(val_all))
  expect_true(is.data.frame(val_se))
})


test_that("validate works with a list of input datasets", {
  val_miniTCGA <- validate(miniTCGA, miniRAVmodel)

  expect_true(is.matrix(val_miniTCGA))
  expect_equal(dim(val_miniTCGA), c(17, 4))
})


test_that("Filter out invalid inputs for validate function", {
  expect_error(validate(microTCGA, miniRAVmodel),
               "Provide a study with at least 8 samples.")
  expect_error(validate(miniTCGA, miniRAVmodel, level = "all"),
               "'level = \"all\"' is not available for a list of datasets.")
})


test_that("different arguments", {
  res <- validate(miniTCGA[[1]], miniRAVmodel, method = "spearman",
                  maxFrom = "avgLoading", level = "all")
  
  expect_equal(dim(res), c(17, 8))
  expect_equal(round(res[1,1], 4), 0.1949)
})
