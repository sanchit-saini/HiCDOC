data("exampleHiCDOCDataSet")
data("exampleHiCDOCDataSetProcessed")

test_chromosome_levels <- function(object, expectedLevels){
    expect_identical(levels(object@distances$chromosome), expectedLevels)
    expect_identical(levels(object@selfInteractionRatios$chromosome), expectedLevels)
    expect_identical(seqlevels(object@compartments), expectedLevels)
    expect_identical(seqlevels(object@concordances), expectedLevels)
    expect_identical(seqlevels(object@differences), expectedLevels)
    expect_identical(levels(object@centroids$chromosome), expectedLevels)
}

test_condition_levels <- function(object, expectedLevels){
    expect_identical(levels(object@distances$condition), expectedLevels)
    expect_identical(levels(object@selfInteractionRatios$condition),
                     expectedLevels)
    expect_identical(levels(object@compartments$condition), expectedLevels)
    expect_identical(levels(object@concordances$condition), expectedLevels)
    expect_identical(levels(object@centroids$condition), expectedLevels)
}

test_replicate_levels <- function(object, expectedLevels){
    expect_identical(levels(object@distances$replicate), expectedLevels)
    expect_identical(levels(object@selfInteractionRatios$replicate),
                     expectedLevels)
    expect_identical(levels(object@concordances$replicate), expectedLevels)
}

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
        objectRed <- reduceHiCDOCDataSet(exampleHiCDOCDataSetProcessed, 
                                         chromosomes = "X"),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, "X")
    expect_identical(objectRed@totalBins, c("X" = 120))
    expect_identical(objectRed@weakBins, list("X" = c(171, 200)))
    # Doesn't remove replicates & conditions
    expect_identical(objectRed$replicate,
                     c("R2", "R1", "R1", "R2", "R2", "R1", "R3"))
    expect_identical(objectRed$condition, c("2", "1", "2", "1", "3", "3", "1"))
    # Interactions
    expect_equal(dim(SummarizedExperiment::assay(objectRed)), c(7021, 7))
    expect_equal(round(sum(SummarizedExperiment::assay(objectRed), na.rm=TRUE),2), 
                 35080.62)
    
    # Objects produced by detectCompartments
    test_chromosome_levels(objectRed, "X")
    test_condition_levels(objectRed, c("1", "2", "3"))
    test_replicate_levels(objectRed, c("R1", "R2", "R3"))
})

test_that("reduceHiCDOCDataSet works if select chromosome, keep levels", {
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            exampleHiCDOCDataSetProcessed,
            chromosomes = "X",
            dropLevels = FALSE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, "X")
    expect_identical(objectRed@totalBins, c("X" = 120))
    expect_identical(objectRed@weakBins, list("X" = c(171, 200)))
    # Doesn't remove replicates & conditions
    expect_identical(objectRed$replicate, exampleHiCDOCDataSetProcessed$replicate)
    expect_identical(objectRed$condition, exampleHiCDOCDataSetProcessed$condition)

    # Interactions
    expect_equal(dim(SummarizedExperiment::assay(objectRed)), c(7021, 7))
    expect_equal(round(sum(SummarizedExperiment::assay(objectRed), na.rm=TRUE),2), 
                 35080.62)
    
    # Objects prduced by detectCompartments
    test_chromosome_levels(objectRed, c("X", "Y", "Z"))
    test_condition_levels(objectRed, c("1", "2", "3"))
    test_replicate_levels(objectRed, c("R1", "R2", "R3"))
})

test_that("reduceHiCDOCDataSet works if select condition, drop levels", {
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            exampleHiCDOCDataSetProcessed,
            conditions = c(1, 2),
            dropLevels = TRUE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, exampleHiCDOCDataSetProcessed@chromosomes)
    expect_identical(objectRed@totalBins, exampleHiCDOCDataSetProcessed@totalBins)
    expect_identical(objectRed@weakBins,  exampleHiCDOCDataSetProcessed@weakBins)
    expect_identical(objectRed$replicate, c("R2", "R1", "R1", "R2", "R3"))
    expect_identical(objectRed$condition, c("2", "1", "2", "1", "1"))

    # Interactions
    expect_equal(dim(SummarizedExperiment::assay(objectRed)), c(39524, 5))
    expect_equal(sum(SummarizedExperiment::assay(objectRed), na.rm=TRUE), 
                 122986.4, tolerance=1e-2)
    
    # Objects prduced by detectCompartments
    test_chromosome_levels(objectRed, c("X", "Y", "Z"))
    test_condition_levels(objectRed, c("1", "2"))
    test_replicate_levels(objectRed, c("R1", "R2", "R3"))
})

test_that("reduceHiCDOCDataSet works if select replicate, drop levels", {
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            exampleHiCDOCDataSetProcessed,
            replicate = "R1",
            dropLevels = TRUE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, exampleHiCDOCDataSetProcessed@chromosomes)
    expect_identical(objectRed@totalBins, exampleHiCDOCDataSetProcessed@totalBins)
    expect_identical(objectRed@weakBins,  exampleHiCDOCDataSetProcessed@weakBins)
    expect_identical(objectRed$replicate, c("R1", "R1", "R1"))
    expect_identical(objectRed$condition, c("1", "2", "3"))

    # Interactions
    expect_equal(dim(SummarizedExperiment::assay(objectRed)), c(39524, 3))
    expect_equal(sum(SummarizedExperiment::assay(objectRed), na.rm=TRUE), 
                 51133.36, tolerance=1e-2)
    
    # Objects prduced by detectCompartments
    test_chromosome_levels(objectRed, c("X", "Y", "Z"))
    test_condition_levels(objectRed, c("1", "2", "3"))
    test_replicate_levels(objectRed, c("R1"))
})

test_that("reduceHiCDOCDataSet works if select chr, cond & rep, keep levels", {
    expect_warning(
        objectRed <- reduceHiCDOCDataSet(
            exampleHiCDOCDataSetProcessed,
            chromosome = "X",
            replicate = "R1",
            condition = "1",
            dropLevels = FALSE
        ),
        "You should not reduce a HiCDOCDataSet after"
    )
    # Chromosomes
    expect_identical(objectRed@chromosomes, "X")
    expect_identical(objectRed@totalBins, c("X" = 120))
    expect_identical(objectRed@weakBins,  list("X" = c(171, 200)))
    expect_identical(objectRed$replicate, "R1")
    expect_identical(objectRed$condition, "1")

    # Interactions
    expect_equal(dim(SummarizedExperiment::assay(objectRed)), c(7021, 1))
    expect_equal(sum(SummarizedExperiment::assay(objectRed), na.rm=TRUE), 
                 6998.319, tolerance=1e-2)

    # Objects prduced by detectCompartments
    test_chromosome_levels(objectRed, c("X", "Y", "Z"))
    test_condition_levels(objectRed, c("1", "2", "3"))
    test_replicate_levels(objectRed, c("R1", "R2", "R3"))
})

