test_that("filterSmallChromosomes behave as expected with default filter", {
    object <- HiCDOCExample()
    # No filter on the example dataset
    expect_message(
        object <- filterSmallChromosomes(object),
        "Keeping chromosomes with at least 100 positions."
    )
    expect_equal(length(object@chromosomes), 2)
    expect_equal(nrow(object@interactions), 86736)
})

test_that("filterSmallChromosomes behave as expected with custom filter", {
    object <- HiCDOCExample()
    # Filter on 1 chromosome
    expect_message(
        object <- filterSmallChromosomes(object, 125),
        "Keeping chromosomes with at least 125 positions."
    )
    expect_equal(length(object@chromosomes), 1)
    expect_identical(object@chromosomes, "17")
    expect_equal(nrow(object@interactions), 48768)
    expect_identical(object@parameters$smallChromosomeThreshold, 115)
})
