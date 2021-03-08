test_that("reduceHiCDOCDataSet return correct errors", {
    object <- HiCDOCDataSetExample()
    # On chromosomes
    expect_error(
        reduceHiCDOCDataSet(object, chromosomes = c(3, 4)),
        "Unknown chromosomes"
    )
    expect_error(
        reduceHiCDOCDataSet(object, chromosomes = "chr1"),
        "Unknown chromosome"
    )
    # On conditions
    expect_error(
        reduceHiCDOCDataSet(object, conditions = c(2, 3)),
        "Unknown conditions"
    )
    expect_error(
        reduceHiCDOCDataSet(object, conditions = "cond1"),
        "Unknown condition"
    )
    # On replicates
    expect_error(
        reduceHiCDOCDataSet(object, replicates = c(3, 4)),
        "Unknown replicates"
    )
    expect_error(
        reduceHiCDOCDataSet(object, replicates = "rep1"),
        "Unknown replicate"
    )
})

test_that("reduceHiCDOCDataSet works if select chromosome, dropLevels", {
    object <- HiCDOCDataSetExample()
    # Run a detectCompartments
    object <- detectCompartments(object)
    expect_warning(
        object <- reduceHiCDOCDataSet(object, chromosomes = c("17")),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(object@chromosomes, "17")
    expect_identical(object@totalBins, c("17" = 127))
    expect_identical(object@weakBins, list("17" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(object@replicates, rep(c("1", "2", "3"), 2))
    expect_identical(object@conditions, rep(c("1", "2"), each = 3))
    # Interactions
    expect_identical(levels(object@interactions$chromosome), "17")
    expect_identical(levels(object@interactions$condition), c("1", "2"))
    expect_identical(levels(object@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(object@interactions), 48768)
    # Objects produced by detectCompartments
    expect_identical(levels(object@distances$chromosome), "17")
    expect_identical(levels(object@selfInteractionRatios$chromosome), "17")
    expect_identical(levels(object@compartments$chromosome), "17")
    expect_identical(levels(object@concordances$chromosome), "17")
    expect_identical(levels(object@differences$chromosome), "17")
    expect_identical(levels(object@centroids$chromosome), "17")
    expect_identical(levels(object@positions$chromosome), "17")

    expect_identical(levels(object@distances$condition), c("1", "2"))
    expect_identical(levels(object@selfInteractionRatios$condition), c("1", "2"))
    expect_identical(levels(object@compartments$condition), c("1", "2"))
    expect_identical(levels(object@concordances$condition), c("1", "2"))
    expect_identical(levels(object@centroids$condition), c("1", "2"))

    expect_identical(levels(object@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(object@selfInteractionRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(object@concordances$replicate), c("1", "2", "3"))
})

test_that("reduceHiCDOCDataSet works if select chromosome, keep levels", {
    object <- HiCDOCDataSetExample()
    # Run a detectCompartments
    object <- detectCompartments(object)
    expect_warning(
        object <- reduceHiCDOCDataSet(
            object,
            chromosomes = c("17"),
            dropLevels = FALSE
        ),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(object@chromosomes, "17")
    expect_identical(object@totalBins, c("17" = 127))
    expect_identical(object@weakBins, list("17" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(object@replicates, rep(c("1", "2", "3"), 2))
    expect_identical(object@conditions, rep(c("1", "2"), each = 3))

    # Interactions
    expect_identical(levels(object@interactions$chromosome), c("17", "18"))
    expect_identical(levels(object@interactions$condition), c("1", "2"))
    expect_identical(levels(object@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(object@interactions), 48768)

    # Objects prduced by detectCompartments
    expect_identical(levels(object@distances$chromosome), c("17", "18"))
    expect_identical(levels(object@selfInteractionRatios$chromosome), c("17", "18"))
    expect_identical(levels(object@compartments$chromosome), c("17", "18"))
    expect_identical(levels(object@concordances$chromosome), c("17", "18"))
    expect_identical(levels(object@differences$chromosome), c("17", "18"))
    expect_identical(levels(object@centroids$chromosome), c("17", "18"))
    expect_identical(levels(object@positions$chromosome), c("17", "18"))

    expect_identical(levels(object@distances$condition), c("1", "2"))
    expect_identical(levels(object@selfInteractionRatios$condition), c("1", "2"))
    expect_identical(levels(object@compartments$condition), c("1", "2"))
    expect_identical(levels(object@concordances$condition), c("1", "2"))
    expect_identical(levels(object@centroids$condition), c("1", "2"))

    expect_identical(levels(object@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(object@selfInteractionRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(object@concordances$replicate), c("1", "2", "3"))
})

test_that("reduceHiCDOCDataSet works if select conditions, dropLevels", {
    object <- HiCDOCDataSetExample()
    # Run a detectCompartments
    object <- detectCompartments(object)
    expect_warning(
        object <- reduceHiCDOCDataSet(object, conditions = c("1")),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(object@chromosomes, c("17", "18"))
    expect_identical(object@totalBins, c("17" = 127, "18" = 112))
    expect_identical(object@weakBins, list("17" = NULL, "18" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(object@replicates, c("1", "2", "3"))
    expect_identical(object@conditions, rep(c("1"), each = 3))
    # Interactions
    expect_identical(levels(object@interactions$chromosome), c("17", "18"))
    expect_identical(levels(object@interactions$condition), c("1"))
    expect_identical(levels(object@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(object@interactions), 43368)
    # Objects prduced by detectCompartments
    expect_identical(levels(object@distances$chromosome), c("17", "18"))
    expect_identical(levels(object@selfInteractionRatios$chromosome), c("17", "18"))
    expect_identical(levels(object@compartments$chromosome), c("17", "18"))
    expect_identical(levels(object@concordances$chromosome), c("17", "18"))
    expect_identical(levels(object@differences$chromosome), c("17", "18"))
    expect_identical(levels(object@centroids$chromosome), c("17", "18"))
    expect_identical(levels(object@positions$chromosome), c("17", "18"))

    expect_identical(levels(object@distances$condition), c("1"))
    expect_identical(levels(object@selfInteractionRatios$condition), c("1"))
    expect_identical(levels(object@compartments$condition), c("1"))
    expect_identical(levels(object@concordances$condition), c("1"))
    expect_identical(levels(object@centroids$condition), c("1"))

    expect_identical(levels(object@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(object@selfInteractionRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(object@concordances$replicate), c("1", "2", "3"))
})


test_that("reduceHiCDOCDataSet works if select replicate, dropLevels", {
    object <- HiCDOCDataSetExample()
    # Run a detectCompartments
    object <- detectCompartments(object)
    expect_warning(
        object <- reduceHiCDOCDataSet(object, replicate = c("1")),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(object@chromosomes, c("17", "18"))
    expect_identical(object@totalBins, c("17" = 127, "18" = 112))
    expect_identical(object@weakBins, list("17" = NULL, "18" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(object@replicates, c("1", "1"))
    expect_identical(object@conditions, c("1", "2"))
    # Interactions
    expect_identical(levels(object@interactions$chromosome), c("17", "18"))
    expect_identical(levels(object@interactions$condition), c("1", "2"))
    expect_identical(levels(object@interactions$replicate), c("1"))
    expect_equal(nrow(object@interactions), 28912)
    # Objects prduced by detectCompartments
    expect_identical(levels(object@distances$chromosome), c("17", "18"))
    expect_identical(levels(object@selfInteractionRatios$chromosome), c("17", "18"))
    expect_identical(levels(object@compartments$chromosome), c("17", "18"))
    expect_identical(levels(object@concordances$chromosome), c("17", "18"))
    expect_identical(levels(object@differences$chromosome), c("17", "18"))
    expect_identical(levels(object@centroids$chromosome), c("17", "18"))
    expect_identical(levels(object@positions$chromosome), c("17", "18"))

    expect_identical(levels(object@distances$condition), c("1", "2"))
    expect_identical(levels(object@selfInteractionRatios$condition), c("1", "2"))
    expect_identical(levels(object@compartments$condition), c("1", "2"))
    expect_identical(levels(object@concordances$condition), c("1", "2"))
    expect_identical(levels(object@centroids$condition), c("1", "2"))

    expect_identical(levels(object@distances$replicate), c("1"))
    expect_identical(levels(object@selfInteractionRatios$replicate), c("1"))
    expect_identical(levels(object@concordances$replicate), c("1"))
})

test_that("reduceHiCDOCDataSet works if select chr, cond & rep, keep levels", {
    object <- HiCDOCDataSetExample()
    # Run a detectCompartments
    object <- detectCompartments(object)
    # This case is used in detectCompartments, in parallel mode
    expect_warning(
        object <- reduceHiCDOCDataSet(object,
            chromosomes = "18",
            replicate = "3",
            condition = "2",
            dropLevels = FALSE
        ),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(object@chromosomes, "18")
    expect_identical(object@totalBins, c("18" = 112))
    expect_identical(object@weakBins, list("18" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(object@replicates, c("3"))
    expect_identical(object@conditions, c("2"))
    # Interactions
    expect_identical(levels(object@interactions$chromosome), c("17", "18"))
    expect_identical(levels(object@interactions$condition), c("1", "2"))
    expect_identical(levels(object@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(object@interactions), 6328)
    # Objects prduced by detectCompartments
    expect_identical(levels(object@distances$chromosome), c("17", "18"))
    expect_identical(levels(object@selfInteractionRatios$chromosome), c("17", "18"))
    expect_identical(levels(object@compartments$chromosome), c("17", "18"))
    expect_identical(levels(object@concordances$chromosome), c("17", "18"))
    expect_identical(levels(object@differences$chromosome), c("17", "18"))
    expect_identical(levels(object@centroids$chromosome), c("17", "18"))
    expect_identical(levels(object@positions$chromosome), c("17", "18"))

    expect_identical(levels(object@distances$condition), c("1", "2"))
    expect_identical(levels(object@selfInteractionRatios$condition), c("1", "2"))
    expect_identical(levels(object@compartments$condition), c("1", "2"))
    expect_identical(levels(object@concordances$condition), c("1", "2"))
    expect_identical(levels(object@centroids$condition), c("1", "2"))

    expect_identical(levels(object@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(object@selfInteractionRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(object@concordances$replicate), c("1", "2", "3"))
})
