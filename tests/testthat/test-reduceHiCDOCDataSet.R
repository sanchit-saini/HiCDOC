data(exampleHiCDOCDataSet)
object <- filterSmallChromosomes(exampleHiCDOCDataSet)
object <- filterWeakPositions(object)
object <- filterSparseReplicates(object)
object <- detectCompartments(object)

test_that("reduceHiCDOCDataSet return correct errors", {
    # On chromosomes
    expect_error(
        reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = c(5, 6)),
        "Unknown chromosomes"
    )
    expect_error(
        reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes = "chr1"),
        "Unknown chromosome"
    )
    # On conditions
    expect_error(
        reduceHiCDOCDataSet(exampleHiCDOCDataSet, conditions = c(3, 4)),
        "Unknown condition: 4"
    )
    expect_error(
        reduceHiCDOCDataSet(exampleHiCDOCDataSet, conditions = "cond1"),
        "Unknown condition"
    )
    # On replicates
    expect_error(
        reduceHiCDOCDataSet(exampleHiCDOCDataSet, replicates = c(3, 4)),
        "Unknown replicates"
    )
    expect_error(
        reduceHiCDOCDataSet(exampleHiCDOCDataSet, replicates = "rep1"),
        "Unknown replicate"
    )
})

test_that("reduceHiCDOCDataSet works if select chromosome, dropLevels", {
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(object, chromosomes = c("X")),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, "X")
    expect_identical(objectRed@totalBins, c("X" = 120))
    expect_identical(objectRed@weakBins, list("X" = c(32L, 91L, 120L)))
    # Doesn't remove replicates & conditions
    expect_identical(objectRed@replicates,
                     c("R2", "R1", "R1", "R2", "R2", "R1", "R3"))
    expect_identical(objectRed@conditions, c("2", "1", "2", "1", "3", "3", "1"))
    # Interactions
    expect_identical(levels(objectRed@interactions$chromosome), "X")
    expect_identical(levels(objectRed@interactions$condition), c("1", "2", "3"))
    expect_identical(levels(objectRed@interactions$replicate), c("R1", "R2", "R3"))
    expect_equal(nrow(objectRed@interactions), 34515)
    # Objects produced by detectCompartments
    expect_identical(levels(objectRed@distances$chromosome), "X")
    expect_identical(levels(objectRed@selfInteractionRatios$chromosome), "X")
    expect_identical(levels(objectRed@compartments$chromosome), "X")
    expect_identical(levels(objectRed@concordances$chromosome), "X")
    expect_identical(levels(objectRed@differences$chromosome), "X")
    expect_identical(levels(objectRed@centroids$chromosome), "X")
    expect_identical(levels(objectRed@positions$chromosome), "X")

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
    expect_identical(levels(objectRed@concordances$replicate),
                     expectedLevels)
})

test_that("reduceHiCDOCDataSet works if select chromosome, keep levels", {
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
    expect_identical(objectRed@weakBins, list("X" = c(32L, 91L, 120L)))
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
    expect_equal(nrow(objectRed@interactions), 34515)

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
    object <- filterSparseReplicates(exampleHiCDOCDataSet)
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
                     levels(object@interactions$condition))
    expect_identical(levels(objectRed@interactions$replicate),
                     c("R1"))
    expect_equal(nrow(objectRed@interactions), 51015)

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
    expect_identical(objectRed@weakBins,  list("X" = c(32L, 91L, 120L)))
    expect_identical(objectRed@replicates, "R1")
    expect_identical(objectRed@conditions, "1")

    # Interactions
    expect_identical(levels(objectRed@interactions$chromosome),
                     levels(object@interactions$chromosome))
    expect_identical(levels(objectRed@interactions$condition),
                     levels(object@interactions$condition))
    expect_identical(levels(objectRed@interactions$replicate),
                     levels(object@interactions$replicate))
    expect_equal(nrow(objectRed@interactions), 6903)

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
