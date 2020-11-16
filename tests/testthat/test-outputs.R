context("GenomicSuperSignature functions tests")

data(miniPCAmodel)
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, miniPCAmodel)

test_that("validate.R", {
  expect_true(is(val_all, "data.frame"))
})

test_that("validatedSignatures.R", {
  val_sig <- validatedSignatures(val_all, num.out = 3, scoreCutoff = 0)
  val_sig2 <- validatedSignatures(val_all, num.out = 3, scoreCutoff = 0, indexOnly = TRUE)

  expect_true(is(val_sig, "matrix"))
  expect_true(is(val_sig2, "vector"))
  expect_true(is(val_sig2, "numeric"))
  expect_equal(val_sig2, c(1076, 2538, 338))
})


expect_error(getModel("C3"), "Prior you entered isn\'t available yet.")
expect_error(meshTable(miniPCAmodel, 1))
expect_equal(dim(meshTable(miniPCAmodel, 1994)), c(31, 2))

test_that("findsignature.R", {
  res <- findSignature(miniPCAmodel, keyword = "Bcell")
  res2 <- findSignature(miniPCAmodel, keyword = "Bcell", k = 5)

  expect_true(is(res, "data.frame"))
  expect_true(is(res2, "vector"))
  expect_true(is(res2, "numeric"))
  expect_equal(dim(res), c(2, 2))
  expect_equal(res2, c(16, 17))
})
