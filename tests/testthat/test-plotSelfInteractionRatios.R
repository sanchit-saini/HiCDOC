test_that("plotSelfInteractionRatios returns an error if no compartments", {
    data(exampleHiCDOCDataSet)
    expect_error(
        pp <- plotSelfInteractionRatios(exampleHiCDOCDataSet),
        "No compartments found."
    )
})

test_that("plotSelfInteractionRatios behaves as expected", {
    data(exampleHiCDOCDataSetProcessed)
    expect_error(
        plotSelfInteractionRatios(exampleHiCDOCDataSetProcessed), 
        '"chromosome"')
    expect_error(
        plotSelfInteractionRatios(exampleHiCDOCDataSetProcessed, 4), 
        "Unknown")

    pp <- plotSelfInteractionRatios(exampleHiCDOCDataSetProcessed, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        pp$labels,
        list(
            "colour" = "Compartment",
            "x" = "Compartment",
            "y" = "Interaction difference",
            "title" = paste0("Differences between self-interactions ",
                             "and other interactions"),
            "subtitle" = "Chromosome X"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomPoint")
    expect_is(pp$layers[[2]]$geom, "GeomBoxplot")
    # No error when printed
    expect_error(print(pp), NA)
})
