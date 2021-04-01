test_that("HiCDOCDataSetFromPro works as expected", {
    # Path to each matrix file
    matrixPaths <-
        system.file("extdata", "liver_18_10M_500000.matrix", package="HiCDOC")

    # Path to each bed file
    bedPaths <-
        system.file("extdata", "liver_18_10M_500000.bed", package="HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- c(1)
    conditions <- c(1)

    # Instantiation of data set
    expect_message(
        object <- HiCDOCDataSetFromHiCPro(
            matrixPaths = matrixPaths,
            bedPaths = bedPaths,
            replicates = replicates,
            conditions = conditions
        ),
        "liver_18_10M_500000.matrix")

    expect_equal(nrow(object@interactions), 210)
    expect_identical(object@chromosomes, "18")
    expect_identical(object@conditions, 1)
    expect_identical(object@replicates, 1)
    expect_identical(object@resolution, 500000L)
})
