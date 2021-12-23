data(exampleHiCDOCDataSet)
object <- reduceHiCDOCDataSet(exampleHiCDOCDataSet, chromosomes=c("X", "Z"), conditions = c(1, 2))

test_that("detectCompartments behaves as expected", {
    # Detect Compartments
    set.seed(3215) # Test with 123 : no significant differences
    expect_message(
        object <- detectCompartments(object, parallel = FALSE),
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
    expect_equal(nrow(object@differences), 1)
    expect_equal(
        nrow(object@differences[pvalue.adjusted <= 0.05]),
        0
    )
    expect_is(object@differences$chromosome, "factor")
    expect_is(object@differences$condition.1, "factor")
    expect_is(object@differences$condition.2, "factor")

    # Centroids
    expect_equal(nrow(object@centroids), 6)
    expect_equal(
        sapply(object@centroids$centroid, length),
        c(rep(120, 4), rep(200, 2))
    )
    expect_equal(
        sapply(object@centroids$centroid, mean),
        c(572.776, 631.262, 514.537, 571.618, 676.45, 803.454),
        tolerance = 1e-04
    )
    expect_is(object@centroids$chromosome, "factor")
    expect_is(object@centroids$condition, "factor")
    expect_is(object@centroids$compartment, "factor")
    expect_is(object@centroids$centroid, "list")

    # Compartments
    expect_equal(nrow(object@compartments), 440)
    expect_equal(
        nrow(object@compartments[compartment == "B"]),
        271
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
    expect_equal(nrow(object@concordances), 880)
    expect_equal(
      nrow(object@concordances[compartment == "A"]),
      337
    )
    expect_equal(
      100 * mean(object@concordances$concordance),
      -0.0218919,
      tolerance = 1e-05
    )

    # Distances
    expect_is(object@distances$chromosome, "factor")
    expect_is(object@distances$index, "numeric")
    expect_is(object@distances$condition, "factor")
    expect_is(object@distances$replicate, "factor")
    expect_is(object@distances$compartment, "factor")
    expect_is(object@distances$distance, "numeric")
    expect_equal(nrow(object@distances), 1760)
    expect_equal(mean(object@distances$distance), 5713.632, tolerance = 1e-04)

    # SelfInteractionRatios
    expect_equal(nrow(object@selfInteractionRatios), 879)
    expect_equal(
      mean(object@selfInteractionRatios$ratio),
      350.6047,
      tolerance = 1e-04
    )
    expect_is(object@selfInteractionRatios$chromosome, "factor")
    expect_is(object@selfInteractionRatios$index, "numeric")
    expect_is(object@selfInteractionRatios$condition, "factor")
    expect_is(object@selfInteractionRatios$replicate, "factor")
    expect_is(object@selfInteractionRatios$ratio, "numeric")
})

test_that("detectCompartments behaves as expected in parallel", {
    skip_on_cran("skip parallel tests")
    skip_on_bioc("skip parallel tests")
  # Parallel settings according to OS
  if(.Platform$OS.type == "windows"){
    multiParam <- BiocParallel::SnowParam(workers = 2)
  } else {
    multiParam <- BiocParallel::MulticoreParam(workers = 2)
  }
  BiocParallel::register(multiParam, default = TRUE)

  # Detect Compartments
  expect_message(
    object <- detectCompartments(object, parallel = TRUE),
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
