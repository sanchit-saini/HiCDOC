data(HiCDOCDataSetExample)
object <- filterSparseReplicates(HiCDOCDataSetExample)
object <- filterWeakPositions(object)
object <- reduceHiCDOCDataSet(object, replicates = c("R1", "R2"), conditions = c("1", "2"))

test_that("plotCompartmentChanges behaves as expected", {
    expect_error(
        pp <- plotCompartmentChanges(object),
        "No compartments found."
    )
    set.seed(3215)
    object <- detectCompartments(object, parallel=FALSE)
    expect_error(plotCompartmentChanges(object),
        "l'argument \"chromosome\" est manquant, avec aucune valeur par défaut")
    expect_error(plotCompartmentChanges(object, 5), "Unknown")

    pp <- plotCompartmentChanges(object, 1)
    expect_is(pp, "ggplot")
    expect_is(pp$labels, "list")
    expect_equal(length(pp$labels), 0)
    expect_is(pp$layers[[1]]$geom, "GeomDrawGrob")
    # No error when printed
    expect_error(print(pp), NA)
})
