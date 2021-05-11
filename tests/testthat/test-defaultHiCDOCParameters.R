test_that("defaultHiCDOCParameters has the expected format", {
    expect_is(defaultHiCDOCParameters, "list")
    expect_identical(
        names(defaultHiCDOCParameters),
        c(
            "smallChromosomeThreshold",
            "sparseReplicateThreshold",
            "weakPositionThreshold",
            "loessSampleSize",
            "kMeansDelta",
            "kMeansIterations",
            "kMeansRestarts"
        )
    )
})

test_that("defaultHiCDOCParameters has the expected values", {
    expect_equal(defaultHiCDOCParameters$smallChromosomeThreshold, 100)
    expect_equal(defaultHiCDOCParameters$weakPositionThreshold, 1)
    expect_equal(defaultHiCDOCParameters$sparseReplicateThreshold, 0.3)
    expect_equal(defaultHiCDOCParameters$loessSampleSize, 20000)
    expect_equal(defaultHiCDOCParameters$kMeansDelta, 1e-04)
    expect_equal(defaultHiCDOCParameters$kMeansIterations, 50)
    expect_equal(defaultHiCDOCParameters$kMeansRestarts, 20)
})
