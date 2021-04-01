test_that("plotDistanceEffect behaves as expected", {
    data(exampleHiCDOCDataSet)
    expect_message(pp <- plotDistanceEffect(exampleHiCDOCDataSet))
    expect_is(pp, "ggExtraPlot")
    # No error when printed
    expect_error(print(pp), NA)
})
