data(miniRAVmodel)
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, miniRAVmodel)


test_that("plotAnnotatedPCA function", {
    res <- plotAnnotatedPCA(exprs(dset), miniRAVmodel, PCnum = c(1,2))
    expect_true(all(class(res) == c("gg", "ggplot", "ggarrange")))
    
    expect_error(plotAnnotatedPCA(exprs(dset), miniRAVmodel, PCnum = c(1,9)))
})


test_that("plotValidate function", {
    res <- plotValidate(val_all)
    res2 <- plotValidate(val_all, minClusterSize = 50)
    res3 <- plotValidate(val_all, interactive = TRUE)
    
    expect_true(all(class(res) == c("gg", "ggplot")))
    expect_true(all(class(res2) == c("gg", "ggplot")))
    expect_true(all(class(res3) == c("plotly", "htmlwidget")))
})


test_that("sampleScoreHeatmap function", {
    score <- calculateScore(dset, miniRAVmodel)
    res <- sampleScoreHeatmap(score, 
                              dataName="bcellViper", 
                              modelName="miniRAVmodel")
    expect_true(class(res) == "Heatmap")
})
