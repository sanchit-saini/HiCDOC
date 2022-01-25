test_that("normalizeDistanceEffect behaves as expected", {
    data(exampleHiCDOCDataSet)
    object <- reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = c("X"))
    object <- filterSparseReplicates(object)
    object <- filterWeakPositions(object)
    
    # Apply normalization
    set.seed(123)
    expect_message(norm <- normalizeDistanceEffect(object))
    # Filtering 0 values before normalisation
    expect_equal(length(norm), length(object))
    matAssay <- SummarizedExperiment::assay(norm)
    
    expect_equal(sum(!is.na(matAssay)), 35105)
    expect_equal(round(colSums(matAssay, na.rm=T),3), 
                 c(4400.339, 5384.152, 0, 0, 7968.404, 8152.765, 9181.879))
})
