data(exampleHiCDOCDataSet)

test_that("filterSmallChromosomes behave as expected with default filter", {
    # No filter on the example dataset
    expect_message(
        object <- filterSmallChromosomes(exampleHiCDOCDataSet),
        "Keeping chromosomes with at least 100 positions."
    )
    expect_equal(length(object@chromosomes), 3)
    expect_equal(nrow(object), 40240)
    assay <- SummarizedExperiment::assay(object)
    expect_equal(sum(!is.na(assay)), 166150)
    expect_equal(colSums(assay, na.rm=TRUE), 
                 c(7891910, 20681568, 7810636, 28374121, 
                   15745554, 21976428, 21352885))
})

test_that("filterSmallChromosomes behave as expected with custom filter", {
    # Filter on 1 chromosome
    expect_message(
        object <- filterSmallChromosomes(exampleHiCDOCDataSet, 161),
        "Keeping chromosomes with at least 161 positions."
    )
    expect_identical(object@chromosomes, "Z")
    expect_equal(nrow(object), 20100)
    expect_identical(object@parameters$smallChromosomeThreshold, 161)
    assay <- SummarizedExperiment::assay(object)
    expect_equal(sum(!is.na(assay)), 40200)
    expect_equal(colSums(assay, na.rm=TRUE), 
                 c(0, 0, 0, 22706459, 
                   0, 0, 7635309))
})

