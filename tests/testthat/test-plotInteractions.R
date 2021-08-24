test_that("plotInteractions behaves as expected", {
    data(exampleHiCDOCDataSet)
    object <- reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = c("X", "Y"))
    
    expect_error(plotInteractions(object, 3), "Unknown chromosome")
    expect_error(plotInteractions(object), '"chromosome"')
    pp <- plotInteractions(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        nrow(pp$data),
        nrow(
            object@interactions %>%
            dplyr::filter(chromosome == "X" & interaction > 0)
        )
    )
    expect_equal(pp$labels$title, "Chromosome X")
    expect_equal(pp$labels$x, "")
    expect_equal(pp$labels$y, "")
    expect_equal(pp$labels$z, "interaction")
    expect_equal(pp$labels$fill, "interaction")
    # No error when printed
    expect_error(print(pp), NA)
})
