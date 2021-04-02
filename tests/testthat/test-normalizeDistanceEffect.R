data(exampleHiCDOCDataSet)
object <- filterSparseReplicates(exampleHiCDOCDataSet)
object <- filterWeakPositions(object)

test_that("normalizeDistanceEffect behaves as expected", {
    # Apply normalization
    set.seed(123)
    expect_message(object <- normalizeDistanceEffect(object))
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$interaction, "numeric")
    # Filtering 0 values before normalisation
    expect_equal(nrow(object@interactions), 181566)
    expect_equal(mean(object@interactions$bin.1), 51.43243, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 101.6621, tolerance = 1e-4)
    expect_equal(mean(object@interactions$interaction), 1.002695,
                 tolerance = 1e-4)
})
