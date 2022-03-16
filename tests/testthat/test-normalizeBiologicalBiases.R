test_that("normalizeBiologicalBiases behaves as expected", {
    data(exampleHiCDOCDataSet)
    object <- filterSparseReplicates(exampleHiCDOCDataSet)
    object <- filterWeakPositions(object)
    
    # Apply normalization
    expect_message(
        norm <- normalizeBiologicalBiases(object),
        "Chromosome Z: normalizing biological biases."
    )
    expect_equal(length(norm), length(object))
    assay <- SummarizedExperiment::assay(norm)
    expect_equal(sum(!is.na(assay)), 181566)
    expect_equal(round(colSums(assay, na.rm=TRUE),2), 
                 c(179.82, 179.82, 79.29, 220.81, 
                   179.81, 179.81, 280.61))
})
