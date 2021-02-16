test_that("normalizeDistanceEffect behaves as expected", {
    object <- HiCDOCExample()
    # Apply normalization
    set.seed(123)
    expect_warning(object <- normalizeDistanceEffect(object))
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$value, "numeric")
    # Filtering 0 values before normalisation
    expect_equal(nrow(object@interactions), 86734)
    expect_equal(mean(object@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(object@interactions$value), 0.9997782, tolerance = 1e-4)
})
