test_that("plotAB behaves as expected", {
    object <- HiCDOCExample()
    expect_error(
        pp <- plotAB(object),
        "Please run 'detectCompartments' first."
    )
    set.seed(3215)
    object <- detectCompartments(object)
    expect_error(plotAB(object), 'chromosomeId" est manquant')
    expect_error(plotAB(object, 3), "Unknown")

    pp <- plotAB(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(pp$labels, list(
        "colour" = "Compartment",
        "x" = "Compartment",
        "y" = "Difference of int.",
        "title" = "Chromosome: 17"
    ))
    expect_is(pp$layers[[1]]$geom, "GeomPoint")
    expect_is(pp$layers[[2]]$geom, "GeomBoxplot")
    # No error when printed
    expect_error(print(pp), NA)
})
