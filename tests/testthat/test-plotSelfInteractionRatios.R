data(HiCDOCDataSetExample)
object <- reduceHiCDOCDataSet(HiCDOCDataSetExample, chromosomes = c("X", "Y"))
object <- filterSparseReplicates(object)
object <- filterWeakPositions(object)

test_that("plotSelfInteractionRatios behaves as expected", {
    expect_error(
        pp <- plotSelfInteractionRatios(object),
        "No compartments found."
    )
    set.seed(3215)
    object <- detectCompartments(object)
    expect_error(plotSelfInteractionRatios(object), '"chromosome"')
    expect_error(plotSelfInteractionRatios(object, 3), "Unknown")

    pp <- plotSelfInteractionRatios(object, 1)
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
