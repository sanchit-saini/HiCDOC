object <- HiCDOCDataSetExample()
object <- filterSparseReplicates(object)
object <- filterWeakPositions(object)

test_that("plotCentroids behaves as expected", {
    expect_error(
        pp <- plotCentroids(object),
        "No compartments found."
    )
    set.seed(3215)
    object <- detectCompartments(object, parallel = FALSE)
    expect_error(plotCentroids(object),
        "l'argument \"chromosome\" est manquant, avec aucune valeur par défaut")
    expect_error(plotCentroids(object, 5), "Unknown")

    pp <- plotCentroids(object, 1)
    expect_is(pp, "ggplot")
    expect_identical(
        pp$labels,
        list(
            "x" = "PC1  69.24 %",
            "y" = "PC2  29.92 %",
            "title" = "PCA on centroids of chromosome W",
            "colour" = "group",
            "shape" = "group"
        )
    )
    expect_is(pp$layers[[1]]$geom, "GeomPoint")
    # No error when printed
    expect_error(print(pp), NA)
})
