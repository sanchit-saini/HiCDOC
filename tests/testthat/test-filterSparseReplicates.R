data(exampleHiCDOCDataSet)
test_that("filterSparseReplicates behaves as expected", {
    # default filtering on the example dataset
    expect_message(
        object <- filterSparseReplicates(exampleHiCDOCDataSet),
        "Removed 1 replicate in total",
        fixed = TRUE
    )
    expect_equal(nrow(object), 43480)
    expect_identical(object@chromosomes, c("W", "X", "Y", "Z"))
    expect_identical(object@parameters$sparseReplicateThreshold, 0.3)
    expect_is(object@validAssay, "list")
    expect_equal(object@validAssay, 
                     list("W" = c(1, 2, 4, 5, 6, 7),
                          "X" = c(1, 2, 5, 6, 7),
                          "Y" = 1:7,
                          "Z" = c(4,7)))
    expect_equal(
        mean(SummarizedExperiment::assay(object), na.rm=T),
        729.7549, tolerance = 1e-4)
    expect_equal(
        sum(is.na(SummarizedExperiment::assay(object))),
        118978)
})

test_that("filterSparseReplicates behaves as expected with custom param", {
    # Filter 1 chromosome on the example dataset
    expect_message(
        object <- filterSparseReplicates(exampleHiCDOCDataSet, 0.9995),
        "Removed 4 replicates in total.",
        fixed = TRUE
    )
    expect_equal(nrow(object), 43480)
    expect_identical(object@chromosomes, c("W", "X", "Y", "Z"))
    # expect_identical(object@weakBins, list("18" = NULL))
    expect_identical(object@parameters$sparseReplicateThreshold, 0.9995)
    expect_equal(
        mean(SummarizedExperiment::assay(object), na.rm=T),
        717.1387, tolerance = 1e-4)
    expect_equal(
        sum(is.na(SummarizedExperiment::assay(object))),
        151280)
})
