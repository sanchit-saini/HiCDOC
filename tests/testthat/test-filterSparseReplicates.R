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
    assay <- SummarizedExperiment::assay(object)
    expect_equal(sum(!is.na(assay)), 185382)
    expect_equal(colSums(assay, na.rm=TRUE), 
                 c(9217842, 23119136 , 7810636, 29607663, 
                   17856784, 24182877, 23488484))
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
    assay <- SummarizedExperiment::assay(object)
    expect_equal(sum(!is.na(assay)), 153080)
    expect_equal(colSums(assay, na.rm=TRUE), 
                 c(5338026, 23119136 , 7810636, 29607663, 
                   10799881, 9615769, 23488484))
})
