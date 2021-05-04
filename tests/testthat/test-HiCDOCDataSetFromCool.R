test_that("HiCDOCDataSetFromCool works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M_500000.cool", package = "HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- 1
    conditions <- 1

    # Instantiation of data set
    expect_message(
        object <- HiCDOCDataSetFromCool(
            paths,
            replicates = replicates,
            conditions = conditions
        ),
        "liver_18_10M_500000.cool'")
    expect_equal(nrow(object@interactions), 210)
    expect_identical(object@chromosomes, "18")
    expect_identical(object@conditions, 1)
    expect_identical(object@replicates, 1)
    expect_identical(object@binSize, 500000L)
})


test_that("HiCDOCDataSetFromCool works as expected if mcool", {
    paths <-
        system.file("extdata", "liver_18_10M.mcool", package = "HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- 1
    conditions <- 1

    # Resolution to select in .mcool files
    binSize <- 500000

    # Instantiation of data set
    expect_message(
        object <- HiCDOCDataSetFromCool(
            paths,
            replicates = replicates,
            conditions = conditions,
            binSize = binSize
        ),
        "liver_18_10M.mcool")

    expect_equal(nrow(object@interactions), 210)
    expect_identical(object@chromosomes, "18")
    expect_identical(object@conditions, 1)
    expect_identical(object@replicates, 1)
    expect_identical(object@binSize, 500000L)
})
