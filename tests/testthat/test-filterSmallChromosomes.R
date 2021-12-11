data(exampleHiCDOCDataSet)

test_that("filterSmallChromosomes behave as expected with default filter", {
    # No filter on the example dataset
    expect_message(
        object <- filterSmallChromosomes(exampleHiCDOCDataSet),
        "Keeping chromosomes with at least 100 positions."
    )
    expect_equal(length(object@chromosomes), 3)
    expect_equal(nrow(object), 40240)
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
})
