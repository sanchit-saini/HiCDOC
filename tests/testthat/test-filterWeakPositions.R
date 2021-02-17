test_that("filterWeakPositions behaves as expected", {
    object <- HiCDOCExample()

    # Filter the 0 values in interactions
    expect_message(
        object <- filterWeakPositions(object),
        "Removed 0 positions"
    )
    expect_equal(nrow(object@interactions), 86734)
    expect_identical(object@weakBins, list("17" = NULL, "18" = NULL))

    # Filter values with 50 threshold
    expect_message(
        object <- filterWeakPositions(object, threshold = 50),
        "Removed 3 positions"
    )
    expect_equal(nrow(object@interactions), 84466)
    expect_identical(
        object@weakBins,
        list("17" = c(124L, 125L, 126L), "18" = NULL)
    )
    expect_identical(object@parameters$weakPosThreshold, 50)
})
