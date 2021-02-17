test_that("HiCDOCDefaultParameters have the expected format", {
    expect_is(HiCDOCDefaultParameters, "list")
    expect_identical(
        names(HiCDOCDefaultParameters),
        c(
            "minLengthChr", "weakPosThreshold", "sparseThreshold",
            "sampleSize", "kMeansIterations", "kMeansDelta",
            "kMeansRestarts"
        )
    )
})

test_that("HiCDOCDefaultParameters have the expected values", {
    expect_equal(HiCDOCDefaultParameters$minLengthChr, 100)
    expect_equal(HiCDOCDefaultParameters$weakPosThreshold, 0)
    expect_equal(HiCDOCDefaultParameters$sparseThreshold, 0.95)
    expect_equal(HiCDOCDefaultParameters$sampleSize, 20000)
    expect_equal(HiCDOCDefaultParameters$kMeansIterations, 50)
    expect_equal(HiCDOCDefaultParameters$kMeansDelta, 1e-04)
    expect_equal(HiCDOCDefaultParameters$kMeansRestarts, 20)
})
