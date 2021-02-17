object <- HiCDOCExample()

test_that("normalizeTechnicalBiases behaves as expected", {
    # Apply normalization
    expect_warning(obj <- normalizeTechnicalBiases(object))
    # Keep obj format
    expect_is(obj@interactions$chromosome, "factor")
    expect_is(obj@interactions$bin.1, "integer")
    expect_is(obj@interactions$bin.2, "integer")
    expect_is(obj@interactions$condition, "factor")
    expect_is(obj@interactions$replicate, "factor")
    expect_is(obj@interactions$value, "numeric")
    # Keep 0 values before normalisation
    expect_equal(nrow(obj@interactions), 86736)
    expect_equal(mean(obj@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$value), 303.4614, tolerance = 1e-4)
})


test_that("normalizeTechnicalBiases behaves as expected in parallel", {
    multiParam <- BiocParallel::MulticoreParam(workers = 3)
    BiocParallel::register(multiParam, default = TRUE)
    # Apply normalization
    obj <- normalizeTechnicalBiases(object, parallel = TRUE)
    # Keep obj format
    expect_is(obj@interactions$chromosome, "factor")
    expect_is(obj@interactions$bin.1, "integer")
    expect_is(obj@interactions$bin.2, "integer")
    expect_is(obj@interactions$condition, "factor")
    expect_is(obj@interactions$replicate, "factor")
    expect_is(obj@interactions$value, "numeric")
    # Keep 0 values before normalisation
    expect_equal(nrow(obj@interactions), 86736)
    expect_equal(mean(obj@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(obj@interactions$value), 303.4614, tolerance = 1e-4)
})
