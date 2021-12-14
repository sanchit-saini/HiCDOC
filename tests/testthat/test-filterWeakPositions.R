data(exampleHiCDOCDataSet)
test_that("filterWeakPositions behaves as expected", {
    # Filter the 0 values in interactions
    expect_message(
        object <- filterWeakPositions(exampleHiCDOCDataSet),
        "Removed 6 positions in total."
    )
    expect_equal(nrow(object), 42447)
    expect_identical(object@weakBins,
                     list("W" = numeric(0),
                          "X" = c(112, 171, 200),
                          "Y" = c(201, 322, 346),
                          "Z" = numeric(0)))
})

test_that("filterWeakPositions behaves as expected, custom param", {
    # Filter values with 50 threshold
    expect_message(
        object <- filterWeakPositions(exampleHiCDOCDataSet, threshold = 50),
        "Removed 123 positions in total."
    )
    expect_equal(nrow(object), 26301)
    expect_identical(sapply(object@weakBins, length),
                     c("W" = 0L, "X" = 120L, "Y" = 3L, "Z" = 0L))
    expect_identical(object@parameters$weakPositionThreshold, 50)
})
