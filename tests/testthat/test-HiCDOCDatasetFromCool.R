test_that("HiCDOCDataSetFromCool produce correct format", {
    path <- system.file(
        "extdata",
        "liver_18_10M_500000.cool",
        package = "HiCDOC"
    )
    replicates <- c("R1", "R2", "R1", "R1", "R2")
    conditions <- c(1, 1, 2, 3, 3)
    expect_error(
        object <- HiCDOCDataSetFromCool(
            rep(path, 5),
            replicates = replicates,
            conditions = conditions,
        ),
        NA
    )
    # Class and slots
    expect_is(object, "HiCDOCDataSet")
    expect_identical(
        slotNames(object),
        c("input", "parameters", "interactions", "chromosomes", "conditions",
          "replicates", "positions", "binSize", "totalBins", "weakBins",
          "validConditions", "validReplicates", "compartments",
          "concordances", "differences", "distances", "centroids",
          "selfInteractionRatios"
        )
    )
    # Class of slots
    expect_is(object@input, "character")
    expect_is(object@interactions, "tbl_df")
    # expect_is(object@weakBins, "list")
    expect_is(object@chromosomes, "character")
    expect_is(object@replicates, "character")
    expect_is(object@conditions, "numeric")
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
    expect_true(is.numeric(object@interactions$interaction))
    # Positions
    expect_is(object@positions$chromosome, "factor")
    expect_is(object@positions$bin, "integer")
    expect_true(is.numeric(object@positions$start))
    expect_true(is.numeric(object@positions$end))
})
