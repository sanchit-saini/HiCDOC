test_that("plotConcordanceDifferences returns an error of no compartments", {
    data(exampleHiCDOCDataSet)
    expect_error(
        pp <- plotConcordanceDifferences(exampleHiCDOCDataSet),
        "Missing slots: comparisons"
    )
})

test_that("plotConcordanceDifferences behaves as expected", {
    data(exampleHiCDOCDataSetProcessed)
    expect_error(
        pp <- plotConcordanceDifferences(exampleHiCDOCDataSetProcessed),
        NA
    )
    expect_is(pp, "ggplot")
    expect_identical(
        unlist(pp$labels),
        c(
            "x" = "Concordance",
            "title" = "Distribution of concordance differences",
            "fill" = "changed",
            "y" = "count",
            "weight" = "weight"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomBar")
    # No error when printed
    expect_error(print(pp), NA)
})
