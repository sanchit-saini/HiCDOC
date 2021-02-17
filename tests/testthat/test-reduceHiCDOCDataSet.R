test_that("reduceHiCDOCDataSet return correct errors", {
    object <- HiCDOCExample()
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
    obj <- HiCDOCExample()
    # Run a detectCompartments
    obj <- detectCompartments(obj)
    expect_warning(
        obj <- reduceHiCDOCDataSet(obj,
            chromosomes = c("17")
        ),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(obj@chromosomes, "17")
    expect_identical(obj@totalBins, c("17" = 127))
    expect_identical(obj@weakBins, list("17" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(obj@replicates, rep(c("1", "2", "3"), 2))
    expect_identical(obj@conditions, rep(c("1", "2"), each = 3))
    # Interactions
    expect_identical(levels(obj@interactions$chromosome), "17")
    expect_identical(levels(obj@interactions$condition), c("1", "2"))
    expect_identical(levels(obj@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(obj@interactions), 48768)
    # Objects produced by detectCompartments
    expect_identical(levels(obj@distances$chromosome), "17")
    expect_identical(levels(obj@diagonalRatios$chromosome), "17")
    expect_identical(levels(obj@compartments$chromosome), "17")
    expect_identical(levels(obj@concordances$chromosome), "17")
    expect_identical(levels(obj@differences$chromosome), "17")
    expect_identical(levels(obj@centroids$chromosome), "17")
    expect_identical(levels(obj@positions$chromosome), "17")

    expect_identical(levels(obj@distances$condition), c("1", "2"))
    expect_identical(levels(obj@diagonalRatios$condition), c("1", "2"))
    expect_identical(levels(obj@compartments$condition), c("1", "2"))
    expect_identical(levels(obj@concordances$condition), c("1", "2"))
    expect_identical(levels(obj@centroids$condition), c("1", "2"))

    expect_identical(levels(obj@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@diagonalRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@concordances$replicate), c("1", "2", "3"))
})

test_that("reduceHiCDOCDataSet works if select chromosome, keep levels", {
    obj <- HiCDOCExample()
    # Run a detectCompartments
    obj <- detectCompartments(obj)
    expect_warning(
        obj <- reduceHiCDOCDataSet(obj,
            chromosomes = c("17"),
            dropLevels = FALSE
        ),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(obj@chromosomes, "17")
    expect_identical(obj@totalBins, c("17" = 127))
    expect_identical(obj@weakBins, list("17" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(obj@replicates, rep(c("1", "2", "3"), 2))
    expect_identical(obj@conditions, rep(c("1", "2"), each = 3))

    # Interactions
    expect_identical(levels(obj@interactions$chromosome), c("17", "18"))
    expect_identical(levels(obj@interactions$condition), c("1", "2"))
    expect_identical(levels(obj@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(obj@interactions), 48768)

    # Objects prduced by detectCompartments
    expect_identical(levels(obj@distances$chromosome), c("17", "18"))
    expect_identical(levels(obj@diagonalRatios$chromosome), c("17", "18"))
    expect_identical(levels(obj@compartments$chromosome), c("17", "18"))
    expect_identical(levels(obj@concordances$chromosome), c("17", "18"))
    expect_identical(levels(obj@differences$chromosome), c("17", "18"))
    expect_identical(levels(obj@centroids$chromosome), c("17", "18"))
    expect_identical(levels(obj@positions$chromosome), c("17", "18"))

    expect_identical(levels(obj@distances$condition), c("1", "2"))
    expect_identical(levels(obj@diagonalRatios$condition), c("1", "2"))
    expect_identical(levels(obj@compartments$condition), c("1", "2"))
    expect_identical(levels(obj@concordances$condition), c("1", "2"))
    expect_identical(levels(obj@centroids$condition), c("1", "2"))

    expect_identical(levels(obj@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@diagonalRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@concordances$replicate), c("1", "2", "3"))
})

test_that("reduceHiCDOCDataSet works if select conditions, dropLevels", {
    obj <- HiCDOCExample()
    # Run a detectCompartments
    obj <- detectCompartments(obj)
    expect_warning(
        obj <- reduceHiCDOCDataSet(obj, conditions = c("1")),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(obj@chromosomes, c("17", "18"))
    expect_identical(obj@totalBins, c("17" = 127, "18" = 112))
    expect_identical(obj@weakBins, list("17" = NULL, "18" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(obj@replicates, c("1", "2", "3"))
    expect_identical(obj@conditions, rep(c("1"), each = 3))
    # Interactions
    expect_identical(levels(obj@interactions$chromosome), c("17", "18"))
    expect_identical(levels(obj@interactions$condition), c("1"))
    expect_identical(levels(obj@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(obj@interactions), 43368)
    # Objects prduced by detectCompartments
    expect_identical(levels(obj@distances$chromosome), c("17", "18"))
    expect_identical(levels(obj@diagonalRatios$chromosome), c("17", "18"))
    expect_identical(levels(obj@compartments$chromosome), c("17", "18"))
    expect_identical(levels(obj@concordances$chromosome), c("17", "18"))
    expect_identical(levels(obj@differences$chromosome), c("17", "18"))
    expect_identical(levels(obj@centroids$chromosome), c("17", "18"))
    expect_identical(levels(obj@positions$chromosome), c("17", "18"))

    expect_identical(levels(obj@distances$condition), c("1"))
    expect_identical(levels(obj@diagonalRatios$condition), c("1"))
    expect_identical(levels(obj@compartments$condition), c("1"))
    expect_identical(levels(obj@concordances$condition), c("1"))
    expect_identical(levels(obj@centroids$condition), c("1"))

    expect_identical(levels(obj@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@diagonalRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@concordances$replicate), c("1", "2", "3"))
})


test_that("reduceHiCDOCDataSet works if select replicate, dropLevels", {
    obj <- HiCDOCExample()
    # Run a detectCompartments
    obj <- detectCompartments(obj)
    expect_warning(
        obj <- reduceHiCDOCDataSet(obj, replicate = c("1")),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(obj@chromosomes, c("17", "18"))
    expect_identical(obj@totalBins, c("17" = 127, "18" = 112))
    expect_identical(obj@weakBins, list("17" = NULL, "18" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(obj@replicates, c("1", "1"))
    expect_identical(obj@conditions, c("1", "2"))
    # Interactions
    expect_identical(levels(obj@interactions$chromosome), c("17", "18"))
    expect_identical(levels(obj@interactions$condition), c("1", "2"))
    expect_identical(levels(obj@interactions$replicate), c("1"))
    expect_equal(nrow(obj@interactions), 28912)
    # Objects prduced by detectCompartments
    expect_identical(levels(obj@distances$chromosome), c("17", "18"))
    expect_identical(levels(obj@diagonalRatios$chromosome), c("17", "18"))
    expect_identical(levels(obj@compartments$chromosome), c("17", "18"))
    expect_identical(levels(obj@concordances$chromosome), c("17", "18"))
    expect_identical(levels(obj@differences$chromosome), c("17", "18"))
    expect_identical(levels(obj@centroids$chromosome), c("17", "18"))
    expect_identical(levels(obj@positions$chromosome), c("17", "18"))

    expect_identical(levels(obj@distances$condition), c("1", "2"))
    expect_identical(levels(obj@diagonalRatios$condition), c("1", "2"))
    expect_identical(levels(obj@compartments$condition), c("1", "2"))
    expect_identical(levels(obj@concordances$condition), c("1", "2"))
    expect_identical(levels(obj@centroids$condition), c("1", "2"))

    expect_identical(levels(obj@distances$replicate), c("1"))
    expect_identical(levels(obj@diagonalRatios$replicate), c("1"))
    expect_identical(levels(obj@concordances$replicate), c("1"))
})

test_that("reduceHiCDOCDataSet works if select chr, cond & rep, keep levels", {
    obj <- HiCDOCExample()
    # Run a detectCompartments
    obj <- detectCompartments(obj)
    # This case is used in detectCompartments, in parallel mode
    expect_warning(
        obj <- reduceHiCDOCDataSet(obj,
            chromosomes = "18",
            replicate = "3",
            condition = "2",
            dropLevels = FALSE
        ),
        "You should not reduce an HiCDOCDataSet object after"
    )
    # Chromosomes
    expect_identical(obj@chromosomes, "18")
    expect_identical(obj@totalBins, c("18" = 112))
    expect_identical(obj@weakBins, list("18" = NULL))
    # Doesn't remove replicates & conditions
    expect_identical(obj@replicates, c("3"))
    expect_identical(obj@conditions, c("2"))
    # Interactions
    expect_identical(levels(obj@interactions$chromosome), c("17", "18"))
    expect_identical(levels(obj@interactions$condition), c("1", "2"))
    expect_identical(levels(obj@interactions$replicate), c("1", "2", "3"))
    expect_equal(nrow(obj@interactions), 6328)
    # Objects prduced by detectCompartments
    expect_identical(levels(obj@distances$chromosome), c("17", "18"))
    expect_identical(levels(obj@diagonalRatios$chromosome), c("17", "18"))
    expect_identical(levels(obj@compartments$chromosome), c("17", "18"))
    expect_identical(levels(obj@concordances$chromosome), c("17", "18"))
    expect_identical(levels(obj@differences$chromosome), c("17", "18"))
    expect_identical(levels(obj@centroids$chromosome), c("17", "18"))
    expect_identical(levels(obj@positions$chromosome), c("17", "18"))

    expect_identical(levels(obj@distances$condition), c("1", "2"))
    expect_identical(levels(obj@diagonalRatios$condition), c("1", "2"))
    expect_identical(levels(obj@compartments$condition), c("1", "2"))
    expect_identical(levels(obj@concordances$condition), c("1", "2"))
    expect_identical(levels(obj@centroids$condition), c("1", "2"))

    expect_identical(levels(obj@distances$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@diagonalRatios$replicate), c("1", "2", "3"))
    expect_identical(levels(obj@concordances$replicate), c("1", "2", "3"))
})
