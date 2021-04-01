test_that("plotDistanceEffect behaves as expected", {
    data(HiCDOCDataSetExample)
    expect_message(pp <- plotDistanceEffect(HiCDOCDataSetExample))
    expect_is(pp, "ggExtraPlot")
    # No error when printed
    expect_error(print(pp), NA)
})
