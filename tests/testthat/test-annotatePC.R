context("Test annotatePC function")

data(miniRAVmodel)
library(bcellViper)
data(bcellViper)

val_all <- validate(dset, miniRAVmodel)


test_that("Users should choose one of the top 8 PCs of their dataset", {
  res <- annotatePC(2, val_all, miniRAVmodel)
  expect_equal(names(res), "PC2.RAV1076")

  res2 <- annotatePC(2, val_all, miniRAVmodel, trimed_pathway_len = 10)
  expect_equal(res2[2,], "KEGG_SPLIC...")
  expect_error(annotatePC(9, val_all, miniRAVmodel))
  expect_error(annotatePC(c(1:2, 9), val_all, miniRAVmodel))
})

test_that("Extract multiple PCs", {
  res3 <- annotatePC(2:4, val_all, miniRAVmodel, simplify = FALSE)
  expect_true(is.list(res3))
  expect_equal(names(res3), c("PC2-RAV1076", "PC3-noAnnot", "PC4-noAnnot"))

  res4 <- annotatePC(2:4, val_all, miniRAVmodel, n = 7)
  expect_true(is.data.frame(res4))
  expect_equal(dim(res4), c(7,3))
})
