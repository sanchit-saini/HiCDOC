test_that("plotInteractions behaves as expected", {
    object <- HiCDOCDataSetExample()
    expect_error(plotInteractions(object, 3), "Unknown chromosome")
    expect_error(plotInteractions(object), '"chromosomeId"')
    pp <- plotInteractions(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        nrow(pp$data),
        nrow(
            object@interactions %>%
            dplyr::filter(chromosome == "17" & value > 0)
        )
    )
    expect_equal(pp$labels$title, "Chromosome: 17")
    expect_equal(pp$labels$x, "")
    expect_equal(pp$labels$y, "")
    # No error when printed
    expect_error(print(pp), NA)
})
