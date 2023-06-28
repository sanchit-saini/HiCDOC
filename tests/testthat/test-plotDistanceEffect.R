test_that("plotDistanceEffect behaves as expected", {
    data(exampleHiCDOCDataSet)
    object <- reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = "X")
    expect_message(pp <- plotDistanceEffect(exampleHiCDOCDataSet))
    expect_is(pp, "ggplot")
    # No error when printed
    expect_error(print(pp), NA)
})
