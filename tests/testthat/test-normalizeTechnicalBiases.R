data(exampleHiCDOCDataSet)
object <- filterSparseReplicates(exampleHiCDOCDataSet)
object <- filterWeakPositions(object)

test_that("normalizeTechnicalBiases behaves as expected", {
    # Apply normalization
    expect_warning(object <- normalizeTechnicalBiases(object, parallel = FALSE))
    # Keep object format
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_is(object@interactions$interaction, "numeric")
    # Keep 0 values before normalisation
    expect_equal(nrow(object@interactions), 181566)
    expect_equal(mean(object@interactions$bin.1), 51.43243, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 101.6621, tolerance = 1e-4)
    expect_equal(mean(object@interactions$interaction), 343.0371,
                 tolerance = 1e-4)
})


test_that("normalizeTechnicalBiases behaves as expected in parallel", {
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
    expect_is(object@interactions$interaction, "numeric")
    # Keep 0 values before normalisation
    expect_equal(nrow(object@interactions), 181566)
    expect_equal(mean(object@interactions$bin.1), 51.43243, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 101.6621, tolerance = 1e-4)
    expect_equal(mean(object@interactions$interaction), 343.0371,
                 tolerance = 1e-4)
})
