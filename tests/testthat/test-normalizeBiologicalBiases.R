obj <- HiCDOCExample()

test_that("normalizeBiologicalBiases behaves as expected", {
    # Apply normalization
    expect_message(
        obj <- normalizeBiologicalBiases(obj),
        "Chromosome: 18"
    )
    # Keep object format
    expect_is(obj@interactions$chromosome, "factor")
    expect_is(obj@interactions$bin.1, "integer")
    expect_is(obj@interactions$bin.2, "integer")
    expect_is(obj@interactions$condition, "factor")
    expect_is(obj@interactions$replicate, "factor")
    expect_is(obj@interactions$value, "numeric")
    # Remove 2 0 values before normalization
    expect_equal(nrow(obj@interactions), 86734)
    expect_equal(mean(obj@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$value), 0.009349513, tolerance = 1e-7)
})
