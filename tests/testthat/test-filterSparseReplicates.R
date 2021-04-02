data(exampleHiCDOCDataSet)
test_that("filterSparseReplicates behaves as expected", {
    # default filtering on the example dataset
    expect_message(
        object <- filterSparseReplicates(exampleHiCDOCDataSet),
        "Removed 8 replicates in total",
        fixed = TRUE
    )
    expect_equal(nrow(object@interactions), 185382)
    expect_identical(object@chromosomes, c("W", "X", "Y", "Z"))
    expect_identical(object@parameters$sparseReplicateThreshold, 0.05)
    expect_is(object@validConditions, "list")
    expect_identical(object@validConditions$W,
                     c('1', '1', '1', '2', '3', '3'))
    expect_identical(object@validConditions$X,
                     c('1', '1', '2', '3', '3'))
    expect_identical(object@validConditions$Y,
                     c('1', '1', '1', '2', '2', '3', '3'))
    expect_identical(object@validConditions$Z,
                     c('1', '1'))
    expect_is(object@validReplicates, "list")
    expect_identical(object@validReplicates$W,
                     c('R1','R2','R3','R2','R1','R2'))
    expect_identical(object@validReplicates$X,
                     c('R1','R3','R2','R1','R2'))
    expect_identical(object@validReplicates$Y,
                     c('R1','R2','R3','R1','R2','R1','R2'))
    expect_identical(object@validReplicates$Z,
                     c('R2','R3'))
})

test_that("filterSparseReplicates behaves as expected with custom param", {
    # Filter 1 chromosome on the example dataset
    expect_message(
        object <- filterSparseReplicates(exampleHiCDOCDataSet, 0.9995),
        "Removed 11 replicates in total.",
        fixed = TRUE
    )
    expect_equal(nrow(object@interactions), 153080)
    expect_identical(object@chromosomes, c("W", "X", "Y", "Z"))
    # expect_identical(object@weakBins, list("18" = NULL))
    expect_identical(object@parameters$sparseReplicateThreshold, 0.9995)
})
