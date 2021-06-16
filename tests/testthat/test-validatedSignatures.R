data(miniRAVmodel)
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, miniRAVmodel)

test_that("validatedSignatures.R", {
    res <- validatedSignatures(val_all, num.out = 3, scoreCutoff = 0)
    res2 <- validatedSignatures(val_all, num.out = 3, scoreCutoff = 0,
                                indexOnly = TRUE)
    res3 <- validatedSignatures(val_all, num.out = 10, scoreCutoff = 0.575)
    
    expect_true(is.matrix(res))
    expect_true(is.vector(res2))
    expect_true(is.numeric(res2))
    expect_equal(res2, c(1076, 2538, 338))
    expect_equal(dim(res3), c(2, 5))
})
