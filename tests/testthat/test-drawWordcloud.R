data(miniRAVmodel)

test_that("meshTable function", {
    res <- meshTable(miniRAVmodel, 1139)
    expect_equal(res[1,1], "NOTCH1 protein, human")
    expect_equal(round(res[1,2], 4), 0.1242)

    res2 <- meshTable(miniRAVmodel, 1139, weighted = FALSE)
    expect_equal(res2[1,1], "Muromonab-CD3")
    expect_equal(round(res2[1,2], 4), 0.3333)

    res3 <- meshTable(miniRAVmodel, 1139, weighted = FALSE, rm.noise = 5)
    expect_equal(res3[1,1], "Infant")
    expect_equal(round(res3[1,2], 4), 0.1429)

    expect_error(meshTable(miniRAVmodel, 1),
                 "Selected ind \\(RAV1\\) doesn't exist.")
    expect_equal(dim(meshTable(miniRAVmodel, 1994)), c(31, 2))
})

expect_error(drawWordcloud(miniRAVmodel, 1),
             "Selected ind \\(RAV1\\) doesn't exist.")
