context("Test getRAVInfo and getStudyinfo functions")

data(miniRAVmodel)

test_that("getRAVInfo", {
    res <- getRAVInfo(miniRAVmodel, ind = 438)
    val <- c(5,7,9,3,5,6)

    expect_true(is.list(res))
    expect_equal(length(res), 4)
    expect_true(all(unlist(res[1:3]) == c(6, 0.04, 43)))
    expect_equal(res[[4]]$PC, val)
})

test_that("getStudyInfo", {
    res2 <- getStudyInfo(miniRAVmodel, "SRP119465")
    val2 <- c("studyTitle", "studySize", "RAVs")

    expect_equal(names(res2), val2)
    expect_equal(dim(res2$RAVs), c(3, 3))
    expect_equal(res2$studySize, 100)
})
