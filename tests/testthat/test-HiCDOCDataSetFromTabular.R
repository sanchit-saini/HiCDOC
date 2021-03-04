test_that("HiCDOCDataSetFromTabular produce correct format", {
    linkToMatrix <- system.file(
        "extdata",
        "sample.tsv",
        package = "HiCDOC"
    )
    expect_error(object <- HiCDOCDataSetFromTabular(linkToMatrix), NA)
    # Class and slots
    expect_is(object, "HiCDOCDataSet")
    expect_identical(
        slotNames(object),
        c(
            "input",
            "parameters",
            "interactions",
            "chromosomes",
            "replicates",
            "positions",
            "conditions",
            "binSize",
            "totalBins",
            "weakBins",
            "compartments",
            "concordances",
            "differences",
            "distances",
            "centroids",
            "selfInteractionRatios",
        )
    )
    # Class of slots
    expect_is(object@input, "character")
    expect_is(object@interactions, "tbl_df")
    expect_is(object@weakBins, "list")
    expect_is(object@chromosomes, "character")
    expect_is(object@replicates, "character")
    expect_is(object@conditions, "character")
    expect_is(object@totalBins, "numeric")
    expect_is(object@binSize, "integer")
    expect_is(object@distances, "NULL")
    expect_is(object@selfInteractionRatios, "NULL")
    expect_is(object@compartments, "NULL")
    expect_is(object@concordances, "NULL")
    expect_is(object@differences, "NULL")
    expect_is(object@centroids, "NULL")
    expect_is(object@parameters, "list")
    expect_is(object@positions, "tbl_df")
    # Interactions
    expect_is(object@interactions$chromosome, "factor")
    expect_is(object@interactions$bin.1, "integer")
    expect_is(object@interactions$bin.2, "integer")
    expect_is(object@interactions$condition, "factor")
    expect_is(object@interactions$replicate, "factor")
    expect_true(is.numeric(object@interactions$value))
    # Positions
    expect_is(object@positions$chromosome, "factor")
    expect_is(object@positions$bin, "integer")
    expect_true(is.numeric(object@positions$start))
    expect_true(is.numeric(object@positions$end))
})

test_that("HiCDOCDataSetFromTabular produce correct values", {
    linkToMatrix <- system.file(
        "extdata",
        "sample.tsv",
        package = "HiCDOC"
    )
    expect_error(object <- HiCDOCDataSetFromTabular(linkToMatrix), NA)

    # Interactions
    expect_equal(nrow(object@interactions), 86736)
    expect_equal(mean(object@interactions$bin.1), 40.81129, tolerance = 1e-5)
    expect_equal(mean(object@interactions$bin.2), 80.62258, tolerance = 1e-5)
    expect_equal(mean(object@interactions$value), 320.2741, tolerance = 1e-5)
    # weakBins
    expect_identical(names(object@weakBins), c("17", "18"))
    # chromosomes
    expect_identical(object@chromosomes, c("17", "18"))
    # replicates & conditions
    expect_identical(object@replicates, c("1", "2", "3", "1", "2", "3"))
    expect_identical(object@conditions, c("1", "1", "1", "2", "2", "2"))
    # bins
    expect_identical(object@totalBins, c("17" = 127, "18" = 112))
    expect_identical(object@binSize, 500000L)
    # Parameters
    expect_identical(object@parameters, defaultHiCDOCParameters)
    # Positions
    expect_equal(mean(object@positions$bin), 60.48536, tolerance = 1e-5)
    expect_equal(mean(object@positions$start), 29742678)
    expect_equal(mean(object@positions$end), 30242677)
})
