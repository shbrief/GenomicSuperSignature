context("Test validate function handles differnet types of inputs")

data(miniRAVmodel)
library(bcellViper)
data(bcellViper)
miniTCGA <- readRDS("miniTCGA.rds")
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

expect_error(validate(microTCGA, miniRAVmodel),
             "Provide a study with at least 8 samples.")
