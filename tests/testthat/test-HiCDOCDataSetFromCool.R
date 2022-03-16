test_that("HiCDOCDataSetFromCool works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M_500000.cool", package = "HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- "GG"
    conditions <- "AA"

    # Instantiation of data set
    expect_message(
        object <- HiCDOCDataSetFromCool(
            paths,
            replicates = replicates,
            conditions = conditions
        ),
        "liver_18_10M_500000.cool'")
    matAssay <- SummarizedExperiment::assay(object) 
    expect_equal(dim(matAssay), c(210, 1))
    expect_identical(object@chromosomes, "18")
    expect_identical(object$replicate, c("GG"))
    expect_identical(object$condition, c("AA"))
})


test_that("HiCDOCDataSetFromCool works as expected if mcool", {
    paths <-
        system.file("extdata", "liver_18_10M.mcool", package = "HiCDOC")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- "A"
    conditions <- "C"

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
    
    matAssay <- SummarizedExperiment::assay(object) 
    expect_equal(dim(matAssay), c(210, 1))
    expect_identical(object@chromosomes, "18")
    expect_identical(object$replicate, c("A"))
    expect_identical(object$condition, c("C"))
})
