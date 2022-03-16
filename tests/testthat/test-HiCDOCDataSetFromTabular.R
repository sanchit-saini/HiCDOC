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
        c("input", "parameters", "chromosomes", 
          "totalBins", "weakBins", "validAssay",
          "compartments", "concordances", "differences", "comparisons", 
          "distances", "centroids", "selfInteractionRatios", "interactions", 
          "colData", "assays", "NAMES", "elementMetadata", "metadata"
        )
    )
    # Class of slots
    expect_is(object@input, "character")
    expect_is(object@weakBins, "list")
    expect_is(object@validAssay, "list")
    expect_is(object@chromosomes, "character")
    expect_is(object@totalBins, "numeric")
    expect_is(object@distances, "NULL")
    expect_is(object@selfInteractionRatios, "NULL")
    expect_is(object@compartments, "NULL")
    expect_is(object@concordances, "NULL")
    expect_is(object@differences, "NULL")
    expect_is(object@centroids, "NULL")
    expect_is(object@parameters, "list")
    expect_is(object, "InteractionSet")
    
    # Interactions
    expect_is(SummarizedExperiment::assay(object), "matrix")
    expect_is(InteractionSet::regions(object), "GRanges")
    expect_is(InteractionSet::interactions(object), "StrictGInteractions")
    expect_is(S4Vectors::mcols(object), "DataFrame")
    expect_true(is.numeric(SummarizedExperiment::assay(object)))
})

test_that("HiCDOCDalinkToMatrixtaSetFromTabular produce correct values", {
    path <- system.file(
        "extdata",
        "liver_18_10M_500000.tsv",
        package = "HiCDOC"
    )
    expect_error(object <- HiCDOCDataSetFromTabular(path), NA)
    gi <- InteractionSet::interactions(object)
    assays <- SummarizedExperiment::assay(object)
    
    # Interactions
    expect_equal(nrow(assays), 210)
    expect_equal(mean(gi@anchor1), 7.333333, tolerance = 1e-5)
    expect_equal(mean(gi@anchor2), 13.66667, tolerance = 1e-5)
    expect_equal(mean(assays), 484.019, tolerance = 1e-5)
    
    # chromosomes
    expect_identical(object@chromosomes, "18")
    # bins
    expect_equal(object@totalBins, c("18" = 20))
    # Parameters
    expect_identical(object@parameters, defaultHiCDOCParameters)
    # Positions
    regions <- data.frame(InteractionSet::regions(object))
    expect_equal(mean(regions$index), 10.5, tolerance = 1e-5)
    expect_equal(mean(regions$start), 4750000)
    expect_equal(mean(regions$end), 5249999)
})
