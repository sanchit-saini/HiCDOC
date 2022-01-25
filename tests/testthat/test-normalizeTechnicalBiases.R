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
                 c(751294.607, 1152226.832, 0, 0, 721966.582, 721539.165, 910899.793))
})
