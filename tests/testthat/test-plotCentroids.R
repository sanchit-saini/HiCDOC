test_that("plotCentroids returns error if no centroids", {
    data(exampleHiCDOCDataSet)
    expect_error(
        pp <- plotCentroids(exampleHiCDOCDataSet),
        "No compartments found."
    )
})

test_that("plotCentroids behaves as expected", {
    data(exampleHiCDOCDataSetProcessed)
    expect_error(plotCentroids(exampleHiCDOCDataSetProcessed),
        "argument \"chromosome\"")
    expect_error(plotCentroids(exampleHiCDOCDataSetProcessed, 5), "Unknown")

    pp <- plotCentroids(exampleHiCDOCDataSetProcessed, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        unlist(pp$labels),
        c("caption" = "Quality controls:\nCentroid PC1 inertia: OK\nA/B clustering consistency: OK",
          "x" = "PC1  91.27 %",
          "y" = "PC2  6.76 %",
          "title" = "PCA on centroids of chromosome X",
          "colour" = "compartment",
          "shape" = "condition"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomPoint")
    # No error when printed
    expect_error(print(pp), NA)
})
