test_that("normalizeTechnicalBiases behaves as expected", {
    data(exampleHiCDOCDataSet)
    object <- reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = c("X"))
    object <- filterSparseReplicates(object)
    object <- filterWeakPositions(object)
    
    # Apply normalization
    expect_warning(norm <- normalizeTechnicalBiases(object, parallel = FALSE))
    # Keep object format
    expect_equal(nrow(norm), nrow(object))
    assay <- SummarizedExperiment::assay(norm)
    expect_equal(sum(!is.na(assay)), 35105)
    expect_equal(round(colSums(assay, na.rm=TRUE),2), 
                 c(1152264.61, 0, 910932.08, 0, 751282.49, 
                   721601.27, 722004.73))
})
