context("Test subsetEnrichedPathways function")

data(miniRAVmodel)

test_that("Confirm the output structure", {
    res <- subsetEnrichedPathways(miniRAVmodel, ind = 695)
    val <- c("IRIS_Bcell-Memory_IgG_IgA", "DMAP_BCELLA3",
             "IRIS_Bcell-naive", "IRIS_Bcell-Memory_IgM", "DMAP_BCELLA2")

    expect_equal(dim(res), c(10, 1))
    expect_equal(colnames(res), "RAV695.Description")
    expect_equal(res[1:5,], val)

    res2 <- subsetEnrichedPathways(miniRAVmodel, ind = 695, include_nes = TRUE)
    expect_true(is.na(res2[10,1]))
})

test_that("Test `include_nes` argument", {
    res3 <- subsetEnrichedPathways(miniRAVmodel, ind = c(695, 1994),
                                   n = 3, include_nes = TRUE)
    val2 <- c("RAV695.Description", "RAV695.NES",
              "RAV1994.Description", "RAV1994.NES")

    expect_equal(colnames(res3), val2)
    expect_equal(dim(res3), c(3, 4))
    expect_equal(round(res3[,4], digits = 6), c(2.682385, 2.676483, 2.634790))
})
