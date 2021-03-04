test_that("normalizeTechnicalBiases behaves as expected", {
    object <- HiCDOCExample()
    # Apply normalization
    expect_warning(object <- normalizeTechnicalBiases(object, parallel = FALSE))
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$value, "numeric")
    # Keep 0 values before normalisation
    expect_equal(nrow(object@interactions), 86736)
    expect_equal(mean(object@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(object@interactions$value), 303.4614, tolerance = 1e-4)
})


test_that("normalizeTechnicalBiases behaves as expected in parallel", {
    object <- HiCDOCExample()
    multiParam <- BiocParallel::MulticoreParam(workers = 3)
    BiocParallel::register(multiParam, default = TRUE)
    # Apply normalization
    object <- normalizeTechnicalBiases(object, parallel = TRUE)
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$value, "numeric")
    # Keep 0 values before normalisation
    expect_equal(nrow(object@interactions), 86736)
    expect_equal(mean(object@interactions$bin.1), 40.81129, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 80.62258, tolerance = 1e-4)
    expect_equal(mean(object@interactions$value), 303.4614, tolerance = 1e-4)
})
