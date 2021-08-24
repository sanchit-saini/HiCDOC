test_that("plotCompartments returns error if no compartments", {
    data(exampleHiCDOCDataSet)
    expect_error(
        pp <- plotCompartments(exampleHiCDOCDataSet),
        "No compartments found."
    )
})

test_that("plotCompartments behaves as expected", {
    data(exampleHiCDOCDataSetProcessed)
    expect_error(plotCompartments(exampleHiCDOCDataSetProcessed),
        "argument \"chromosome\"")
    expect_error(plotCompartments(exampleHiCDOCDataSetProcessed, 5), "Unknown")

    pp <- plotCompartments(exampleHiCDOCDataSetProcessed, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        unlist(pp$labels),
        c(
          "x" = "position",
          "fill" = "compartment",
          "y" = "count",
          "weight" = "weight"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomBar")
    # No error when printed
    expect_error(print(pp), NA)
})
