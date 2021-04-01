data(HiCDOCDataSetExample)

test_that("reduceHiCDOCDataSet return correct errors", {
    # On chromosomes
    expect_error(
        reduceHiCDOCDataSet(HiCDOCDataSetExample, chromosomes = c(5, 6)),
        "Unknown chromosomes"
    )
    expect_error(
        reduceHiCDOCDataSet(HiCDOCDataSetExample, chromosomes = "chr1"),
        "Unknown chromosome"
    )
    # On conditions
    expect_error(
        reduceHiCDOCDataSet(HiCDOCDataSetExample, conditions = c(3, 4)),
        "Unknown condition: 4"
    )
    expect_error(
        reduceHiCDOCDataSet(HiCDOCDataSetExample, conditions = "cond1"),
        "Unknown condition"
    )
    # On replicates
    expect_error(
        reduceHiCDOCDataSet(HiCDOCDataSetExample, replicates = c(3, 4)),
        "Unknown replicates"
    )
    expect_error(
        reduceHiCDOCDataSet(HiCDOCDataSetExample, replicates = "rep1"),
        "Unknown replicate"
    )
})

test_that("reduceHiCDOCDataSet works if select chromosome, dropLevels", {
    # Run a detectCompartments
    object <- filterSparseReplicates(HiCDOCDataSetExample)
    object <- filterWeakPositions(object)
    object <- detectCompartments(object)
    expect_warning(
        object <- reduceHiCDOCDataSet(object, chromosomes = c("X")),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(object@chromosomes, "X")
    expect_identical(object@totalBins, c("X" = 120))
    expect_identical(object@weakBins, list("X" = c(91L, 120L)))
    # Doesn't remove replicates & conditions
    expect_identical(object@replicates,
                     c("R2", "R1", "R1", "R2", "R2", "R1", "R3"))
    expect_identical(object@conditions, c("2", "1", "2", "1", "3", "3", "1"))
    # Interactions
    expect_identical(levels(object@interactions$chromosome), "X")
    expect_identical(levels(object@interactions$condition), c("1", "2", "3"))
    expect_identical(levels(object@interactions$replicate), c("R1", "R2", "R3"))
    expect_equal(nrow(object@interactions), 35105)
    # Objects produced by detectCompartments
    expect_identical(levels(object@distances$chromosome), "X")
    expect_identical(levels(object@selfInteractionRatios$chromosome), "X")
    expect_identical(levels(object@compartments$chromosome), "X")
    expect_identical(levels(object@concordances$chromosome), "X")
    expect_identical(levels(object@differences$chromosome), "X")
    expect_identical(levels(object@centroids$chromosome), "X")
    expect_identical(levels(object@positions$chromosome), "X")

    expect_identical(levels(object@distances$condition), c("1", "2", "3"))
    expect_identical(levels(object@selfInteractionRatios$condition),
                     c("1", "2", "3"))
    expect_identical(levels(object@compartments$condition), c("1", "2", "3"))
    expect_identical(levels(object@concordances$condition), c("1", "2", "3"))
    expect_identical(levels(object@centroids$condition), c("1", "2", "3"))

    expect_identical(levels(object@distances$replicate), c("R1", "R2", "R3"))
    expect_identical(levels(object@selfInteractionRatios$replicate),
                     c("R1", "R2", "R3"))
    expect_identical(levels(object@concordances$replicate),
                     c("R1", "R2", "R3"))
})

test_that("reduceHiCDOCDataSet works if select chromosome, keep levels", {
    object <- filterSparseReplicates(HiCDOCDataSetExample)
    object <- filterWeakPositions(object)
    # Run a detectCompartments
    object <- detectCompartments(object)
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            object,
            chromosomes = c("X"),
            dropLevels = FALSE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, "X")
    expect_identical(objectRed@totalBins, c("X" = 120))
    expect_identical(objectRed@weakBins, list("X" = c(91L, 120L)))
    # Doesn't remove replicates & conditions
    expect_identical(objectRed@replicates, object@replicates)
    expect_identical(objectRed@conditions, object@conditions)

    # Interactions
    expect_identical(levels(objectRed@interactions$chromosome),
                     levels(object@interactions$chromosome))
    expect_identical(levels(objectRed@interactions$condition),
                     levels(object@interactions$condition))
    expect_identical(levels(objectRed@interactions$replicate),
                     levels(object@interactions$replicate))
    expect_equal(nrow(objectRed@interactions), 35105)

    expectedLevels <- levels(object@interactions$chromosome)
    # Objects prduced by detectCompartments
    expect_identical(levels(objectRed@distances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@differences$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@centroids$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@positions$chromosome),
                     expectedLevels)

    expectedLevels <- levels(object@interactions$condition)
    expect_identical(levels(objectRed@distances$condition), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$condition),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$condition), expectedLevels)
    expect_identical(levels(objectRed@concordances$condition), expectedLevels)
    expect_identical(levels(objectRed@centroids$condition), expectedLevels)

    expectedLevels <- levels(object@interactions$replicate)
    expect_identical(levels(objectRed@distances$replicate), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$replicate),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$replicate), expectedLevels)
})

test_that("reduceHiCDOCDataSet works if select condition, drop levels", {
    object <- filterSparseReplicates(HiCDOCDataSetExample)
    object <- filterWeakPositions(object)
    # Run a detectCompartments
    object <- detectCompartments(object, parallel=FALSE)
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            object,
            conditions = c(1, 2),
            dropLevels = TRUE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, object@chromosomes)
    expect_identical(objectRed@totalBins, object@totalBins)
    expect_identical(objectRed@weakBins,  object@weakBins)
    expect_identical(objectRed@replicates, c("R2", "R1", "R1", "R2", "R3"))
    expect_identical(objectRed@conditions, c("2", "1", "2", "1", "1"))

    # Interactions
    expect_identical(levels(objectRed@interactions$chromosome),
                     levels(object@interactions$chromosome))
    expect_identical(levels(objectRed@interactions$condition),
                     c("1", "2"))
    expect_identical(levels(objectRed@interactions$replicate),
                     levels(object@interactions$replicate))
    expect_equal(nrow(objectRed@interactions), 136238)

    expectedLevels <- levels(object@interactions$chromosome)
    # Objects prduced by detectCompartments
    expect_identical(levels(objectRed@distances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@differences$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@centroids$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@positions$chromosome),
                     expectedLevels)

    expectedLevels <- c("1", "2")
    expect_identical(levels(objectRed@distances$condition), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$condition),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$condition), expectedLevels)
    expect_identical(levels(objectRed@concordances$condition), expectedLevels)
    expect_identical(levels(objectRed@centroids$condition), expectedLevels)

    expectedLevels <- levels(object@interactions$replicate)
    expect_identical(levels(objectRed@distances$replicate), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$replicate),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$replicate), expectedLevels)
})

test_that("reduceHiCDOCDataSet works if select replicate, drop levels", {
    object <- filterSparseReplicates(HiCDOCDataSetExample)
    object <- filterWeakPositions(object)
    # Run a detectCompartments
    object <- detectCompartments(object, parallel=FALSE)
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            object,
            replicate = c("R1"),
            dropLevels = TRUE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, object@chromosomes)
    expect_identical(objectRed@totalBins, object@totalBins)
    expect_identical(objectRed@weakBins,  object@weakBins)
    expect_identical(objectRed@replicates, c("R1", "R1", "R1"))
    expect_identical(objectRed@conditions, c("1", "2", "3"))

    # Interactions
    expect_identical(levels(objectRed@interactions$chromosome),
                     levels(object@interactions$chromosome))
    expect_identical(levels(objectRed@interactions$condition),
                     levels(object@conditions))
    expect_identical(levels(objectRed@interactions$replicate),
                     c("R1"))
    expect_equal(nrow(objectRed@interactions), 57731)

    expectedLevels <- levels(object@interactions$chromosome)
    # Objects prduced by detectCompartments
    expect_identical(levels(objectRed@distances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@differences$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@centroids$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@positions$chromosome),
                     expectedLevels)

    expectedLevels <- levels(object@interactions$condition)
    expect_identical(levels(objectRed@distances$condition), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$condition),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$condition), expectedLevels)
    expect_identical(levels(objectRed@concordances$condition), expectedLevels)
    expect_identical(levels(objectRed@centroids$condition), expectedLevels)

    expectedLevels <- c("R1")
    expect_identical(levels(objectRed@distances$replicate), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$replicate), expectedLevels)
    expect_identical(levels(objectRed@concordances$replicate), expectedLevels)
})

test_that("reduceHiCDOCDataSet works if select chr, cond & rep, keep levels", {
    # case used in detectCompartments
    object <- filterSparseReplicates(HiCDOCDataSetExample)
    object <- filterWeakPositions(object)
    # Run a detectCompartments
    object <- detectCompartments(object, parallel=FALSE)
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            object,
            chromosome = c("X"),
            replicate = c("R1"),
            condition = c("1"),
            dropLevels = FALSE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, "X")
    expect_identical(objectRed@totalBins, c("X" = 120))
    expect_identical(objectRed@weakBins,  list("X" = c(91L, 120L)))
    expect_identical(objectRed@replicates, "R1")
    expect_identical(objectRed@conditions, "1")

    # Interactions
    expect_identical(levels(objectRed@interactions$chromosome),
                     levels(object@interactions$chromosome))
    expect_identical(levels(objectRed@interactions$condition),
                     levels(object@interactions$condition))
    expect_identical(levels(objectRed@interactions$replicate),
                     levels(object@interactions$replicate))
    expect_equal(nrow(objectRed@interactions), 7021)

    expectedLevels <- levels(object@interactions$chromosome)
    # Objects prduced by detectCompartments
    expect_identical(levels(objectRed@distances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@differences$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@centroids$chromosome),
                     expectedLevels)
    expect_identical(levels(objectRed@positions$chromosome),
                     expectedLevels)

    expectedLevels <- levels(object@interactions$condition)
    expect_identical(levels(objectRed@distances$condition), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$condition),
                     expectedLevels)
    expect_identical(levels(objectRed@compartments$condition), expectedLevels)
    expect_identical(levels(objectRed@concordances$condition), expectedLevels)
    expect_identical(levels(objectRed@centroids$condition), expectedLevels)

    expectedLevels <- levels(object@interactions$replicate)
    expect_identical(levels(objectRed@distances$replicate), expectedLevels)
    expect_identical(levels(objectRed@selfInteractionRatios$replicate),
                     expectedLevels)
    expect_identical(levels(objectRed@concordances$replicate), expectedLevels)
})
