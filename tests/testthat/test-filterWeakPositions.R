data(exampleHiCDOCDataSet)
test_that("filterWeakPositions behaves as expected", {
    # Filter the 0 values in interactions
    expect_message(
        object <- filterWeakPositions(exampleHiCDOCDataSet),
        "Removed 6 positions in total."
    )
    expect_equal(nrow(object@interactions), 181180)
    expect_identical(object@weakBins,
                     list("W" = NULL,
                          "X" = c(32L, 91L, 120L),
                          "Y" = c(1L, 122L, 146L),
                          "Z" = NULL))
})

test_that("filterWeakPositions behaves as expected, custom param", {
    # Filter values with 50 threshold
    expect_message(
        object <- filterWeakPositions(exampleHiCDOCDataSet, threshold = 50),
        "Removed 123 positions in total."
    )
    expect_equal(nrow(object@interactions), 146467)
    expect_identical(sapply(object@weakBins, length),
                     c("W" = 0L, "X" = 120L, "Y" = 3L, "Z" = 0L))
    expect_identical(object@parameters$weakPositionThreshold, 50)
})
