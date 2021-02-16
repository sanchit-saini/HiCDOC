test_that("normalizeBiologicalBiases behaves as expected", {
    object <- HiCDOCExample()
    # Apply normalization
    expect_message(object <- normalizeBiologicalBiases(object), 
                   "Chromosome: 18")
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$value, "numeric")
    # Remove 2 0 values before normalization
    expect_equal(nrow(object@interactions), 86734)
    expect_equal(mean(object@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(object@interactions$value), 0.009349513, tolerance = 1e-7)
})
