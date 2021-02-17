test_that("plotDiffConcordances behaves as expected", {
    object <- HiCDOCExample()
    expect_error(
        pp <- plotDiffConcordances(object),
        "Please run 'detectCompartments' first."
    )
    set.seed(3215)
    object <- detectCompartments(object)
    pp <- plotDiffConcordances(object)
    expect_is(pp, "ggplot")
    expect_identical(
        pp$labels,
        list(
            "x" = "Concordance",
            "title" = "Distribution of the differences of concordances",
            "fill" = "changed",
            "y" = "count",
            "weight" = "weight"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomBar")
    # No error when printed
    expect_error(print(pp), NA)
})
