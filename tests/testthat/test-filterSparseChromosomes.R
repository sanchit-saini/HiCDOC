object <- HiCDOCExample()

test_that("filterSparseChromosomes behaves as expected", {
    # No filtering on the example dataset
    expect_message(filterSparseChromosomes(object),
        "No chromosome removed (threshold: 95%)",
        fixed = TRUE
    )
    expect_identical(object@chromosomes, c("17", "18"))
    expect_identical(object@parameters$sparseThreshold, 0.95)

    # Filter 1 chromosome on the example dataset
    expect_message(object <- filterSparseChromosomes(object, 0.0001),
        "1 chromosome(s) removed: 17 (threshold : 0.01%)",
        fixed = TRUE
    )
    expect_equal(nrow(object@interactions), 37968)
    expect_identical(object@chromosomes, c("18"))
    expect_identical(object@weakBins, list("18" = NULL))
    expect_identical(object@parameters$sparseThreshold, 0.0001)
})
