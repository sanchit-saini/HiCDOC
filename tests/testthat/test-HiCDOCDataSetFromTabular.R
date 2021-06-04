test_that("HiCDOCDataSetFromTabular produce correct format", {
    path <- system.file(
        "extdata",
        "liver_18_10M_500000.tsv",
        package = "HiCDOC"
    )
    expect_error(object <- HiCDOCDataSetFromTabular(path), NA)
    # Class and slots
    expect_is(object, "HiCDOCDataSet")
    expect_identical(
        slotNames(object),
        c("input", "parameters", "interactions", "chromosomes", "conditions",
          "replicates", "positions", "binSize", "totalBins", "weakBins",
          "validConditions", "validReplicates", "compartments",
          "concordances", "differences", "comparisons", "distances", "centroids",
          "selfInteractionRatios"
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
    expect_true(is.numeric(object@interactions$interaction))
    # Positions
    expect_is(object@positions$chromosome, "factor")
    expect_is(object@positions$bin, "integer")
    expect_true(is.numeric(object@positions$start))
    expect_true(is.numeric(object@positions$end))
})

test_that("HiCDOCDalinkToMatrixtaSetFromTabular produce correct values", {
    path <- system.file(
        "extdata",
        "liver_18_10M_500000.tsv",
        package = "HiCDOC"
    )
    expect_error(object <- HiCDOCDataSetFromTabular(path), NA)

    # Interactions
    expect_equal(nrow(object@interactions), 210)
    expect_equal(mean(object@interactions$bin.1), 7.333333, tolerance = 1e-5)
    expect_equal(mean(object@interactions$bin.2), 13.66667, tolerance = 1e-5)
    expect_equal(mean(object@interactions$interaction), 484.019,
                 tolerance = 1e-5)
    # chromosomes
    expect_identical(object@chromosomes, "18")
    # replicates & conditions
    expect_identical(object@replicates, "R1")
    expect_identical(object@conditions, "1")
    # bins
    expect_identical(object@totalBins, c("18" = 20))
    expect_identical(object@binSize, 500000L)
    # Parameters
    expect_identical(object@parameters, defaultHiCDOCParameters)
    # Positions
    expect_equal(mean(object@positions$bin), 10.5, tolerance = 1e-5)
    expect_equal(mean(object@positions$start), 4750000)
    expect_equal(mean(object@positions$end), 5249999)
})
