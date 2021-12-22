data(exampleHiCDOCDataSet)

test_that("detectCompartments behaves as expected", {
    # Detect Compartments
    set.seed(3215) # Test with 123 : no significant differences
    expect_message(
        object <- detectCompartments(exampleHiCDOCDataSet, parallel = FALSE),
        "Detecting significant differences."
    )
    # Create new objects in correct format
    expect_is(object@distances, "data.table")
    expect_is(object@selfInteractionRatios, "data.table")
    expect_is(object@compartments, "data.table")
    expect_is(object@concordances, "data.table")
    expect_is(object@differences, "data.table")
    expect_is(object@centroids, "data.table")

    # Differences
    expect_equal(nrow(object@differences), 52)
    expect_equal(
        nrow(object@differences[pvalue.adjusted <= 0.05]),
        0
    )
    expect_is(object@differences$chromosome, "factor")
    expect_is(object@differences$condition.1, "factor")
    expect_is(object@differences$condition.2, "factor")

    # Centroids
    expect_equal(nrow(object@centroids), 20)
    expect_equal(
        sapply(object@centroids$centroid, length),
        c(rep(80, 6), rep(120, 6), rep(160, 6), rep(200, 2))
    )
    expect_equal(
        sapply(object@centroids$centroid, mean),
        c(
            431.85, 357.606, 562.354, 637.759, 582.983, 702.689, 571.618, 
            514.537, 631.262, 572.776, 914.067, 1036.959, 457.758, 446.806, 
            688.807, 688.714, 857.06, 931.403, 676.45, 803.454
        ),
        tolerance = 1e-04
    )
    expect_is(object@centroids$chromosome, "factor")
    expect_is(object@centroids$condition, "factor")
    expect_is(object@centroids$compartment, "factor")
    expect_is(object@centroids$centroid, "list")

    # Compartments
    expect_equal(nrow(object@compartments), 1280)
    expect_equal(
        nrow(object@compartments[compartment == "B"]),
        745
    )
    expect_is(object@compartments$chromosome, "factor")
    expect_is(object@compartments$condition, "factor")
    expect_is(object@compartments$compartment, "factor")
    expect_is(object@compartments$index, "numeric")

    # Concordance
    expect_is(object@concordances$chromosome, "factor")
    expect_is(object@concordances$index, "numeric")
    expect_is(object@concordances$condition, "factor")
    expect_is(object@concordances$replicate, "factor")
    expect_is(object@concordances$compartment, "factor")
    expect_is(object@concordances$concordance, "numeric")
    expect_equal(nrow(object@concordances), 2720)
    expect_equal(
      nrow(object@concordances[compartment == "A"]),
      1162
    )
    expect_equal(
      100 * mean(object@concordances$concordance),
      0.02854235,
      tolerance = 1e-05
    )

    # Distances
    expect_is(object@distances$chromosome, "factor")
    expect_is(object@distances$index, "numeric")
    expect_is(object@distances$condition, "factor")
    expect_is(object@distances$replicate, "factor")
    expect_is(object@distances$compartment, "factor")
    expect_is(object@distances$distance, "numeric")
    expect_equal(nrow(object@distances), 5440)
    expect_equal(mean(object@distances$distance), 4542.623, tolerance = 1e-04)

    # SelfInteractionRatios
    expect_equal(nrow(object@selfInteractionRatios), 2714)
    expect_equal(
      mean(object@selfInteractionRatios$ratio),
      361.8327,
      tolerance = 1e-04
    )
    expect_is(object@selfInteractionRatios$chromosome, "factor")
    expect_is(object@selfInteractionRatios$index, "numeric")
    expect_is(object@selfInteractionRatios$condition, "factor")
    expect_is(object@selfInteractionRatios$replicate, "factor")
    expect_is(object@selfInteractionRatios$ratio, "numeric")
})

test_that("detectCompartments behaves as expected in parallel", {
  # Parallel settings according to OS
  if(.Platform$OS.type == "windows"){
    multiParam <- BiocParallel::SnowParam(workers = 2)
  } else {
    multiParam <- BiocParallel::MulticoreParam(workers = 2)
  }
  BiocParallel::register(multiParam, default = TRUE)

  # Detect Compartments
  expect_message(
    object <- detectCompartments(exampleHiCDOCDataSet, parallel = TRUE),
    "Detecting significant differences."
  )
  
  # Create new objects in correct format
  expect_is(object@distances, "data.table")
  expect_is(object@selfInteractionRatios, "data.table")
  expect_is(object@compartments, "data.table")
  expect_is(object@concordances, "data.table")
  expect_is(object@differences, "data.table")
  expect_is(object@centroids, "data.table")

  # Differences
  expect_is(object@differences$chromosome, "factor")
  expect_is(object@differences$condition.1, "factor")
  expect_is(object@differences$condition.2, "factor")

  # Centroids
  expect_is(object@centroids$chromosome, "factor")
  expect_is(object@centroids$condition, "factor")
  expect_is(object@centroids$compartment, "factor")
  expect_is(object@centroids$centroid, "list")

  # Compartments
  expect_is(object@compartments$chromosome, "factor")
  expect_is(object@compartments$condition, "factor")
  expect_is(object@compartments$compartment, "factor")
  expect_is(object@compartments$index, "numeric")

  # Concordance
  expect_is(object@concordances$chromosome, "factor")
  expect_is(object@concordances$index, "numeric")
  expect_is(object@concordances$condition, "factor")
  expect_is(object@concordances$replicate, "factor")
  expect_is(object@concordances$compartment, "factor")
  expect_is(object@concordances$concordance, "numeric")

  # Distances
  expect_is(object@distances$chromosome, "factor")
  expect_is(object@distances$index, "numeric")
  expect_is(object@distances$condition, "factor")
  expect_is(object@distances$replicate, "factor")
  expect_is(object@distances$compartment, "factor")
  expect_is(object@distances$distance, "numeric")

  # SelfInteractionRatios
  expect_is(object@selfInteractionRatios$chromosome, "factor")
  expect_is(object@selfInteractionRatios$index, "numeric")
  expect_is(object@selfInteractionRatios$condition, "factor")
  expect_is(object@selfInteractionRatios$replicate, "factor")
  expect_is(object@selfInteractionRatios$ratio, "numeric")
})
