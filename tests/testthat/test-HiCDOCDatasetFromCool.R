test_that("HiCDOCDataSetFromCool produce correct format", {
    basedir <- system.file("extdata", package = "HiCDOC")
    data <- read.csv(file.path(basedir, "coolData.csv"))
    expect_error(
        object <- HiCDOCDataSetFromCool(
            file.path(basedir, data$FileName),
            data$Replicate,
            data$Condition
        ),
        NA
    )
    # Class and slots
    expect_is(object, "HiCDOCDataSet")
    expect_identical(
        slotNames(object),
        c(
            "inputPath",
            "interactions",
            "weakBins",
            "chromosomes",
            "replicates",
            "totalReplicates",
            "totalReplicatesPerCondition",
            "conditions",
            "totalBins",
            "binSize",
            "distances",
            "diagonalRatios",
            "compartments",
            "concordances",
            "differences",
            "centroids",
            "parameters",
            "positions"
        )
    )
    # Class of slots
    expect_is(object@inputPath, "character")
    expect_is(object@interactions, "tbl_df")
    expect_is(object@weakBins, "list")
    expect_is(object@chromosomes, "character")
    expect_is(object@replicates, "character")
    expect_is(object@totalReplicates, "integer")
    expect_is(object@totalReplicatesPerCondition, "numeric")
    expect_is(object@conditions, "integer")
    expect_is(object@totalBins, "numeric")
    expect_is(object@binSize, "integer")
    expect_is(object@distances, "NULL")
    expect_is(object@diagonalRatios, "NULL")
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
