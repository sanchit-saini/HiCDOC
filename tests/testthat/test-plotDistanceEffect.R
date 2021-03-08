test_that("plotDistanceEffect behaves as expected", {
    object <- HiCDOCDataSetExample()
    expect_message(pp <- plotDistanceEffect(object))
    expect_is(pp, "ggExtraPlot")
    # No error when printed
    expect_error(print(pp), NA)
})
