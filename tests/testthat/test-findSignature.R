context("Test functions in findSignature.R file")

data(miniRAVmodel)
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, miniRAVmodel)

test_that("findSignature function", {
    res <- findSignature(miniRAVmodel, keyword = "Bcell")
    res2 <- findSignature(miniRAVmodel, keyword = "Bcell", k = 5)

    expect_true(is.data.frame(res))
    expect_true(is.vector(res2))
    expect_true(is.numeric(res2))
    expect_equal(dim(res), c(2, 2))
    expect_equal(res2, c(16, 17, 18))
    expect_warning(findSignature(miniRAVmodel, keyword = "Bcell", k = 1),
                   "There is no RAV with 1 keyword-containing, enriched pathways.")
})

test_that("findKeywordInRAV function", {
    res <- findKeywordInRAV(miniRAVmodel, "Bcell", ind = 695)
    expect_equal(res, "1|2|3|4|5|6|8 (out of 9)")

    expect_error(findKeywordInRAV(miniRAVmodel, "Bcell", ind =  1),
                 "Selected ind \\(RAV1\\) doesn't exist.")
    expect_warning(findKeywordInRAV(miniRAVmodel, "hello", ind =  695),
                   "RAV695 doesn't have any pathway with the keyword, hello")
})
