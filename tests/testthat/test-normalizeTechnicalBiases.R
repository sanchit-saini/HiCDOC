test_that("normalizeTechnicalBiases behaves as expected", {
    data(exampleHiCDOCDataSet)
    object <- reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = c("X"))
    object <- filterSparseReplicates(object)
    object <- filterWeakPositions(object)
    
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
    expect_equal(nrow(object@interactions), 35105)
    expect_equal(mean(object@interactions$bin.1), 40.0578, tolerance = 1e-4)
    expect_equal(mean(object@interactions$bin.2), 79.4167, tolerance = 1e-4)
    expect_equal(mean(object@interactions$interaction), 852.3587,
                 tolerance = 1e-4)
})
