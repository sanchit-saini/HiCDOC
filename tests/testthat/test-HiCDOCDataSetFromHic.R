test_that("HiCDOCDataSetFromHic works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M.hic", package="HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- c(1)
    conditions <- c(1)
    binSize <- 500000

    # Instantiation of data set
    expect_message(
        object <- HiCDOCDataSetFromHiC(
            paths,
            replicates = replicates,
            conditions = conditions,
            binSize = binSize
        ),
        "liver_18_10M.hic'")
    expect_equal(nrow(object@interactions), 210)
    expect_identical(object@chromosomes, "18")
    expect_identical(object@conditions, 1)
    expect_identical(object@replicates, 1)
    expect_identical(object@binSize, 500000L)
})
