context("Test findStudiesInCluster function")

data(miniRAVmodel)

test_that("studyTitle argument", {
    res <- findStudiesInCluster(miniRAVmodel, 1076)
    res2 <- findStudiesInCluster(miniRAVmodel, 1076, studyTitle = TRUE)

    expect_true(is.character(res))
    expect_true(is.data.frame(res2))
    expect_equal(dim(res2), c(10,2))
    expect_equal(res[1:3], c("SRP028155", "SRP049340", "SRP058840"))
})

test_that("only valid index", {
    expect_error(findStudiesInCluster(miniRAVmodel, 1),
                 "Selected ind \\(RAV1\\) doesn't exist.")
    expect_error(findStudiesInCluster(miniRAVmodel, c(1, 1076, 100)),
                 "Selected ind \\(RAV1, RAV100\\) doesn't exist.")
})
