test_that("HiCDOCDataSetFromHic works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M.hic", package = "HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- 1
    conditions <- 1
    binSize <- 500000

    # Instantiation of data set
    expect_message(
        object <- HiCDOCDataSetFromHiC(
            paths,
            replicates = replicates,
            conditions = conditions,
            binSize = binSize
        ),
        "liver_18_10M.hic")
    expect_equal(nrow(object), 210)
    expect_identical(object@chromosomes, "18")
    expect_identical(object$condition, "X1")
    expect_identical(object$replicat, "1")
})
