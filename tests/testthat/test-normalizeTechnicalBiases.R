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
    
    # Two different possible values, 
    # 1rst with 
    # BLAS:   openblas/libblas.so.3
    # LAPACK: libopenblasp-r0.2.20.so
    # 2nd with 
    # Matrix products: default
    # BLAS:   atlas/libblas.so.3.10.3
    # LAPACK: atlas/liblapack.so.3.10.
    expect_equal(any(sapply(list(
        c(751294.6, 1152226.8, 0, 0, 721966.6, 721539.2, 910899.8),
        c(751294.6, 1152227.4, 0.0, 0.0, 721967.5, 721538.9, 910899.4)), 
        function(x) identical(x, round(colSums(assay, na.rm=TRUE),1)))), TRUE)
    
})
