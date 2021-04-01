data(HiCDOCDataSetExample)
object <- reduceHiCDOCDataSet(HiCDOCDataSetExample,
                              replicates = c("R1", "R2"),
                              conditions = c("1", "2"))
object <- filterSparseReplicates(object)
object <- filterWeakPositions(object)

test_that("plotConcordances behaves as expected", {
    expect_error(
        pp <- plotConcordances(object),
        "No compartments found."
    )
    set.seed(3215)
    object <- detectCompartments(object, parallel=FALSE)
    expect_error(plotConcordances(object), '"chromosome"')
    expect_error(plotConcordances(object, 6), "Unknown")

    pp <- plotConcordances(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        pp$labels,
        list(
            "caption" = "No change is significant (adjusted p-value <= 5%)",
            "x" = "position",
            "y" = "concordance",
            "colour" = "replicate",
            "yintercept" = "yintercept"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomLine")
    # No error when printed
    expect_error(print(pp), NA)
})
