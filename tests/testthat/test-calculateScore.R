data(miniRAVmodel)
data(miniTCGA)
library(bcellViper)
data(bcellViper)

res <- calculateScore(dset, miniRAVmodel)
res2 <- calculateScore(miniTCGA, miniRAVmodel)
res3 <- calculateScore(miniTCGA, miniRAVmodel, rescale.after = FALSE)


test_that("check both input dataset formats works", {
    expect_true(is.matrix(res))
    expect_true(is.list(res2))
})

test_that("check whether `rescale.after` argument works", {
    expect_equal(round(res2[[1]][1,1], 4), -0.7251)
    expect_equal(round(res3[[1]][1,1], 4), -30.822)
})
