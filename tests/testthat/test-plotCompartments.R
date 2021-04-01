data(HiCDOCDataSetExample)
object <- filterSparseReplicates(HiCDOCDataSetExample)
object <- filterWeakPositions(object)

test_that("plotCompartments behaves as expected", {
    expect_error(
        pp <- plotCompartments(object),
        "No compartments found."
    )
    set.seed(3215)
    object <- detectCompartments(object, parallel = FALSE)
    expect_error(plotCompartments(object),
        "argument \"chromosome\"")
    expect_error(plotCompartments(object, 5), "Unknown")

    pp <- plotCompartments(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        pp$labels,
        list(
            "x" = "position",
            "fill" = "compartment",
            "y" = "count",
            "weight" = "weight"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomBar")
    # No error when printed
    expect_error(print(pp), NA)
})
