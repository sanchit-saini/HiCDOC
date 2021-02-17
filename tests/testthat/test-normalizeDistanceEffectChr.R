test_that("normalizeDistanceEffect behaves as expected", {
    obj <- HiCDOCExample()
    # Apply normalization
    set.seed(123)
    expect_warning(obj <- normalizeDistanceEffect(obj))
    # Keep obj format
    expect_is(obj@interactions$chromosome, "factor")
    expect_is(obj@interactions$bin.1, "integer")
    expect_is(obj@interactions$bin.2, "integer")
    expect_is(obj@interactions$condition, "factor")
    expect_is(obj@interactions$replicate, "factor")
    expect_is(obj@interactions$value, "numeric")
    # Filtering 0 values before normalisation
    expect_equal(nrow(obj@interactions), 86734)
    expect_equal(mean(obj@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$value), 0.9997782, tolerance = 1e-4)
})
