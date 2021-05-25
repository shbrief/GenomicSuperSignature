context("Test annotateRAV function")

data(miniRAVmodel)

res <- annotateRAV(miniRAVmodel, ind = 1076)
res2 <- annotateRAV(miniRAVmodel, ind = 1076, abs = TRUE)

expect_equal(res[1,1], "REACTOME_METABOLISM_OF_PROTEINS")
expect_equal(res2[1,1], "REACTOME_CELL_CYCLE")