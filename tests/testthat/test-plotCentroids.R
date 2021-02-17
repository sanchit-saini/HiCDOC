test_that("plotCentroids behaves as expected", {
    object <- HiCDOCExample()
    expect_error(pp <- plotCentroids(object), 
                 "Please run 'detectCompartments' first.")
    set.seed(3215)
    object <- detectCompartments(object)
    expect_error(plotCentroids(object), 'chromosomeId" est manquant')
    expect_error(plotCentroids(object, 3), 'Unknown')
    
    pp <- plotCentroids(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(pp$labels, list("x" = "PC1  94.47 %",
                                     "y" = "PC2  5.19 %",
                                     "title" = "Centroids of chromosome 17",
                                     "colour" = "group",
                                     "shape" = "group"))
    expect_is(pp$layers[[1]]$geom, "GeomPoint")
    # No error when printed
    expect_error(print(pp), NA)
})
